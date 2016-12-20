<?php

include('functions.php');

$cmd_bin = "/usr/local/bin/R --vanilla -e '%s'";
$expression = 'print(system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()),T))';

$stderr = '';
$stdout = '';
$command = sprintf($cmd_bin, $expression);
echo $command;
        $exit_code = cmd_exec($command, $stderr, $stdout);

echo $exit_code;
echo '<pre>';
print_r($stderr);
print_r($stdout);
echo '</pre>';
?>
