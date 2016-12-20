<?php

session_start();

if (($_REQUEST['dataset'] != 'disease') && ($_REQUEST['dataset'] != 'all')  && ($_REQUEST['dataset'] != 'tissue') && ($_REQUEST['dataset'] != 'process') 
		&& ($_REQUEST['dataset'] != 'disease_tissue') && ($_REQUEST['dataset'] != 'disease_process') && ($_REQUEST['dataset'] != 'tissue_process')){
	echo "ERROR:\nWrong query";
}

$folder = 'uploads/';
if (file_exists($folder . session_id() . '/' . $_REQUEST['dataset'] . '/motifs.txt')) {
	echo 1;
} else {
	echo 0;
}

?>