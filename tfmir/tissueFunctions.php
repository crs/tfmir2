<!--done by Maryam -->
<?php


function getTissueOptions($filename) {

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