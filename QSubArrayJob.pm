package QSubArrayJob;
use strict;
use warnings;
use Data::Dumper;
use Log::Log4perl qw(get_logger :levels);
use File::Temp qw/ tempfile tempdir /;
use File::Basename;
use Storable qw(store_fd);
use Bio::SeqIO;

use lib '../modules/';
use FileUtils;

$SIG{'INT'} = \&SIGINT_handler;
$SIG{'QUIT'} = \&SIGQUIT_handler;
$SIG{'ABRT'} = \&SIGABRT_handler;


#logger setup
Log::Log4perl->easy_init($DEBUG);
my $logger = get_logger("QSubUtility");
$logger->level($DEBUG);

my $DEFAULT_JOBS = 400;
my $TIME_LIMIT_SECONDS = 259200;  #roughly 3 days
my @RUNNING_JOBS;
## class (static) variables
our $DEFAULT_TASKCOUNT = 250;
my $DEFAULT_QUEUE = 'standard.q';

#designated qsub executable location
my $QSUB_exec=FileUtils::get_exec("qsub","../config.yml");

=head1 Example Usage

  my $qjob = QSubArrayJob->new();
  # $qjob->queue_name("debug.q");  # default is standard.q
  $qjob->base_dir(dirname($fasta_in));
  $qjob->job_name($job_name);
  $qjob->shard_fasta_input($fasta_in);

  $qjob->task_command(
                      "$EXEC $OPTIONS \$INPUT_FNAME $fasta_db $KTUP > \$OUTPUT_FNAME"
                      );
  $qjob->write_array_job_script();

  $logger->debug("ClusterFASTA: submitting qsub array job");
  $qjob->submit_array_job();

  $logger->debug("Job Number: " . $qjob->job_number());
  $qjob->wait_jobs();

  $logger->debug("ClusterFASTA: jobs complete, collecting output");
  $qjob->collect_output();

  return $qjob->output_fname();

=cut


sub SIGINT_handler {
    # Delete the running jobs, as the script won't be there to reap them
	$logger->info("^C detected, canceling qsub runs");
	map{ print STDERR `qdel $_`}@RUNNING_JOBS;
    exit(0);

}


sub SIGQUIT_handler {
    # Delete the running jobs, as the script won't be there to reap them
	$logger->info("SIG QUIT detected, canceling qsub runs");
	map{ print STDERR `qdel $_`}@RUNNING_JOBS;
    exit(0);

}


sub SIGABRT_handler {
    # Delete the running jobs, as the script won't be there to reap them
	$logger->info("SIG ABRT detected, canceling qsub runs");
	map{ print STDERR `qdel $_`}@RUNNING_JOBS;
    exit(0);

}

# sub SetDefaultQueue {
#   my ($q_name) = @_;
#   $DEFAULT_QUEUE = $q_name;
# }


sub new {
  my ($class, $in_base_dir) = @_;
  my $self = {};
  bless($self, $class);
  $self->{'_DEFAULT_TASKCOUNT'} = \$DEFAULT_TASKCOUNT;
  $self->{'_DEFAULT_QUEUE'} = \$DEFAULT_QUEUE;

  return $self;
}


## private ??
sub count_fasta_sequences {
  my ($self, $fasta_fname) = @_;

  my $sequence_count=`grep -c '>' $fasta_fname`;
  chomp $sequence_count;
  return $sequence_count;
}

=head1 Methods

=cut


=head2 shard_fasta_input

Given the filename of the full fasta sequences, splits input and creates TaskUnits for each shard.

=cut

sub shard_fasta_input {
    my ($self, $full_fname, $JOBS) = @_;

    if (defined($JOBS)) { $self->task_count($JOBS); }
    else                { $self->task_count($self->default_taskcount()); }

    $logger->logconfess("No file found") unless (-e $full_fname && -f $full_fname);

    my $sequence_count = $self->count_fasta_sequences($full_fname);

    $logger->logconfess("No sequences in Fasta file") unless ($sequence_count > 0);

    if ($sequence_count < $self->task_count()) {
      $self->task_count($sequence_count);
    }

    # my $basedir = $self->base_dir();
    # my $tempdir = tempdir ( "tmpseqdir_XXXXX", DIR => $self->base_dir(), CLEANUP=>1);
    $self->temp_dir( tempdir ( "tmpseqdir_XXXXX", DIR => $self->base_dir(), CLEANUP=>1));

    $logger->logconfess("Could not create tmp directory") unless (-e $self->temp_dir() && -d $self->temp_dir());
    my $tmp=[];

    for my $i (1 .. ($self->task_count())) {
      my $tasknum = sprintf("%d", $i);
      my $tfilename = File::Spec->catfile($self->temp_dir(), "input.$tasknum");
      my $toutfilename = File::Spec->catfile($self->temp_dir(), "output.$tasknum");

      my $temp_fh;
      open($temp_fh, ">$tfilename") or $logger->confess("Can't open input shard for writing");
      push @$tmp, {'input_fname'  => $tfilename,
                   'input_fh'     => $temp_fh,
                   'output_fname' => $toutfilename,
                   'status_prefix'=> $tfilename };  # TaskUnit. SHOULD THIS BE AN OBJECT??
    }

    my $seqIO = Bio::SeqIO->new(-format=>'fasta', -file=>$full_fname);
    my $file_i = 0;  # this is index into our $tmp array
    my $db_length = 0;
    while (my $seqobj = $seqIO->next_seq()) {
      my $this_fh = $tmp->[$file_i]->{'input_fh'};
      $db_length += length($seqobj->seq());
      print $this_fh ">" . $seqobj->display_id ."\n" . $seqobj->seq . "\n";
      $file_i++;
      $file_i = $file_i % $self->task_count();  # wrap around
    }

    foreach my $tu (@$tmp) {
      close($tu->{'input_fh'});
      delete $tu->{'input_fh'};  # this file handle should not be considered valid
    }

    $self->task_units($tmp);

    # $self->task_start(1);
    # $self->task_end($self->task_count());
    # return $db_length);
}


=head2 submit_array_job

Runs qsub with the proper options to start all of the tasks.

=cut

sub submit_array_job {
  my ($self, $TIME, $MEM) = @_;

  my $stderr = File::Spec->catfile( $self->temp_dir(), $self->job_name() . ".err");
  my $stdout = File::Spec->catfile( $self->temp_dir(), $self->job_name() . ".out");

  my $qsub_opts = 
    join(" ", ("-r y",
               "-l mem=$MEM,time=$TIME" ,
               "-e $stderr",
               "-o $stdout",
               "-cwd",
               "-N " . $self->job_name(),
               "-t 1-" . $self->task_count()));
  my $cmd_text = "$QSUB_exec $qsub_opts " . $self->job_shell_fname();
  $logger->debug($cmd_text);
  my $cmd_ret = `$cmd_text`;

  $self->job_number( $self->parse_qsub_output($cmd_ret) );
  # ERROR HANDLING WHERE ARE YOU??
}

=head2 parse_qsub_output

Gets the qsub job number from the command line output using a simple regular expression.

=cut

sub parse_qsub_output {
    
  my ($self, $message) = @_;
  $message =~ /Your job-array (\d+)/;
  my $job_id = $1;
    print $job_id, "\n";
  unless (defined $job_id) { $logger->debug( "Cannot find the job id for your the qsub job\n"); exit; }
  return $job_id;
}


=head2 wait_jobs

Filename based polling of jobs.

  $task_unit->{'status_prefix'}.running: still working
  $task_unit->{'status_prefix'}.complete: done

This method blocks until either all jobs finish or a cumulative wait time of
$TIME_LIMIT_SECONDS is reached.

=cut

sub wait_jobs {
  my ($self) = @_;

  foreach my $tu (@{ $self->task_units() }) { $tu->{'complete'} = 0; }

  my $timer = 2;
  my $total_time = 0;
  #use polynomial back-off, as most jobs should complete quickly
  my $all_finished=0;
  until ($all_finished == 1) {
    my @not_finished = grep{ $_->{'complete'} == 0 } @{ $self->task_units() };
    foreach my $tu (@not_finished){
      $tu->{'complete'} = 1 if (-e $tu->{'status_prefix'} . ".complete" &&
                                -f $tu->{'status_prefix'} . ".complete");
    }

    if ( scalar(@not_finished) == 0 || $total_time >= $TIME_LIMIT_SECONDS) {
      $all_finished = 1;
      last;
    }

    $timer+=2 unless $timer > 60;
    $logger->debug("Checking if jobs are complete..." .
                   scalar(@not_finished) .
                   " not finished, pausing $timer s");
    sleep ($timer);
    $total_time += $timer;
  }
}


=head2 collect_output

Concatenates the output files for each task unit into a single file, and setting the output_fname
instance variable.

Also deletes input shards, their respective output files, and status files.

=cut

sub collect_output {
  my ($self) = @_;
    #print "Entering collect_output\n";
    #print $self->temp_dir(), "\n";
  my $coll_out_fname = File::Spec->catfile($self->temp_dir(), 'output.collect');
    #print  $coll_out_fname, "\n";  
  foreach my $hr (@{ $self->task_units() }) {
    my $task_in_fname = $hr->{input_fname};   # this is a temp file
    print $task_in_fname, "\n";
    my $task_out_fname = $hr->{output_fname};    # this should have been created by the process

    $logger->logconfess("No result found for $task_in_fname, looking for filename $task_out_fname")
      unless (-e $task_out_fname && -f $task_out_fname);

    system("cat $task_out_fname >> $coll_out_fname;");
    unlink($hr->{input_fname});
    unlink($hr->{output_fname});
    # unlink($hr->{stderr});
    # unlink($hr->{stdout});
    unlink($hr->{status_prefix}.".complete");
  }
  $logger->logconfess("Some problem in copying the output file, check permissions")
    unless (-e "$coll_out_fname" && -f "$coll_out_fname");

  $self->output_fname($coll_out_fname);
}


=head2 write_array_job_script

Writes the single-task script that will be passed to qsub, in the proper temp directory.

The script as "hard coded" here was designed to take $SGE_TASK_ID environment variable to control
which piece of the array job to execute.
Proper naming of the input shards is taken care of by shard_fasta_input and codified in each
TaskUnit pseudo object.

=cut

sub write_array_job_script {
  my ($self, $fasta_file) = @_;
  my $script_text;
  my $qsub_jobname = $self->job_name();

  # will add a trailing slash if there is none
  my $tempdir_prefix = File::Spec->catfile($self->temp_dir(), '');

  $script_text = <<STARTSCRIPTBLOCK;
#\$ -S /bin/sh
#\$ -N $qsub_jobname
STARTSCRIPTBLOCK

  $script_text .= q{
if [ -z $SGE_TASK_ID ]
then
  echo "SGE_TASK_ID must be set for ArrayJob." 1>&2
  exit
fi
};

  $script_text .= qq{
INPUT_FNAME=${tempdir_prefix}input.\$SGE_TASK_ID
OUTPUT_FNAME=${tempdir_prefix}output.\$SGE_TASK_ID
};


$script_text .= <<STATUSBLOCK;   
echo $fasta_file > \$INPUT_FNAME.running    
date >> \$INPUT_FNAME.running
hostname >> \$INPUT_FNAME.running
pwd >> \$INPUT_FNAME.running
date >> \$INPUT_FNAME.running    

STATUSBLOCK

  # ok here's where we put in the user's cluster task.

  # Is there any reason to use a "pattern" instead of just a string?
  # All commands should be uniform between tasks.
  # $script_text .= sprintf($self->task_command_pattern(), @{ $self->task_command_variables() } );

  $script_text .= $self->task_command();

  $script_text .= <<ENDBLOCK;

mv \$INPUT_FNAME.running \$INPUT_FNAME.complete
#sync

ENDBLOCK

  $self->job_shell_fname( File::Spec->catfile($self->temp_dir(), "task_shell.sh") );

  open (SCRIPTFH, '>'.$self->job_shell_fname()) or $logger->logconfess("Can't write shell command to file system");
  print SCRIPTFH $script_text;
  close (SCRIPTFH);
  return ;
}



=head2 Accessor Functions

These accessors are written in the "perl oo" way. When called without any argument, the accessor
will return the current value of the instance variable (getter). If an argument is provided, the
accessor will set the instance variable in the object (setter).

=cut

#### Accessors ####

## class variable
sub default_queue {
  my ($self, $in_queue_name) = @_;

  if (defined($in_queue_name)) {
    ${ $self->{'_DEFAULT_QUEUE'} } = $in_queue_name;
  }
  return ${ $self->{'_DEFAULT_QUEUE'} };
}


## class variable
sub default_taskcount {
  my ($self, $in_taskcount) = @_;

  if (defined($in_taskcount)) {
    ${ $self->{'_DEFAULT_TASKCOUNT'} } = $in_taskcount;
  }
  return ${ $self->{'_DEFAULT_TASKCOUNT'} };
}

## instance variable
sub queue_name {
  my ($self, $in_qname) = @_;
  if (defined($in_qname)) {
    $self->{'_queue_name'} = $in_qname;
  }

  if (not defined($self->{'_queue_name'})) {
    return $self->default_queue();
  }
  else {
    return $self->{'_queue_name'};
  }
}

## instance variables
sub temp_dir { return hash_accessor('_temp_dir', @_); }
sub job_name { return hash_accessor('_job_name', @_); }
sub base_dir { return hash_accessor('_base_dir', @_); }
sub task_count { return hash_accessor('_task_count', @_); }
sub task_command { return hash_accessor('_task_command', @_); }
sub job_shell_fname { return hash_accessor('_job_shell_fname', @_); }
sub job_number { return hash_accessor('_job_number', @_); }
sub output_fname { return hash_accessor('_output_fname', @_); }

# task_units is a reference to an array of hashes whose specification
# were are loosely referring to as "Task Unit"
# input_fname
# output_fname
# status_prefix  (will mostly be $input_fname, so that we can look for $status_prefix.running/.completed)
sub task_units { return hash_accessor('_task_units', @_); }


# hash_accessor
#
# a "private" function.
# This process was factored out of all of the accessors because of
# the way we are just storing instance variables in a $self hash.
#
sub hash_accessor {
  my ($key, $self, $in_val) = @_;

  if (defined($in_val)) {
    $self->{$key} = $in_val;
  }
  return $self->{$key};
}

1;
