<!--done by Maryam -->
<?php


function getTissueOptions($tissue) {

    
  if ($tissue == "Mus") {
    $filename = "../backend/tissue_Mus.txt";
  }
  
  if ($tissue == "HM") {
    $filename = "../backend/tissue_HM.txt";
  }
  

	$f = fopen($filename, 'r');
	$pb = '';
	if ($f) {
		while (($line = fgets($f)) !== false) {
			$line = trim($line);

			if (mb_detect_encoding($line, 'UTF-8', true)) {
				$pb .= '<option '. checkOption('tissue', $line).' value="'.$line.'">' . $line . '</option>';
			}
		}
	}
	return $pb;
}

?>
<!--done by Maryam -->