<?php


function getDiseaseOptions($species) {
  
  if ($species == "Mus") {
    $filename = "../backend/disease_Mus.txt";
  }
  
  if ($species == "HM") {
    $filename = "../backend/disease_HM.txt";
  }
  
  $f = fopen($filename, 'r');
  $pb = '';
  if ($f) {
    while (($line = fgets($f)) !== false) {
    	$line = trim($line);

    	if (mb_detect_encoding($line, 'UTF-8', true)) {
        	$pb .= '<option '. checkOption('disease', $line).' value="'.$line.'">' . $line . '</option>';
    	}
    }
  }
  return $pb;
}

?>