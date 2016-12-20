<?php header('Content-Type: text/javascript; charset=UTF-8'); 
session_start();
 $url = '//'.$_SERVER['HTTP_HOST'].dirname($_SERVER['PHP_SELF']);

?>


function showMCDS(file) {
	$('#MCDSloader').css('visibility','hidden');
	$("#dialog").dialog({
    			autoOpen: false,
			    position: 'center' ,
			    title: 'MCDS Nodes',
			    draggable: false,
			    resizable : false,
			    modal : true,
			});

	var conDomSetFile = './uploads/<?php echo session_id() ?>/' + file + '/MCDS.txt'
	var loc  = '<?php echo $url; ?>/exists.php?fileName='+conDomSetFile;
	console.log(loc);
	$.get(loc, function(data) {
		if (data == 1) {
			console.log(data);
			//$("#dialog code").load(conDomSetFile, function() {
				//$("#dialog").dialog("open");
			//	console.log(conDomSetFile);
			//});

			$.get(conDomSetFile, function(data) {
				var str = data.split("\n");
				console.log(str);
				var cy = $('.cytoPanel').cytoscape('get');

				cy.elements().unselect();
				for (var i = 0; i < str.length; i++) {

					if (str[i].length == 0) continue;
					console.log(str[i]);
					cy.nodes('[id="' + str[i] + '"]').select();
					cy.nodes('[id="' + str[i].toLowerCase() + '"]').select();
				}
			});
			$('#MCDSloader').css('visibility','hidden');
		} else {
			$('#MCDSloader').css('visibility','visible');
			$.get('<?php echo $url; ?>/cgi-bin/MCDS_scc.cgi?session=<?php echo session_id(); ?>&interactions=' + file, function(conDomSetcompute) {
				console.log("Data not found, try to compute");
				console.log(conDomSetcompute);
				if (conDomSetcompute == 1) {
					console.log("successfully computed connected dominating set");
					showMCDS(file);
				} else {
					alert('Could not compute connected dominating set');
				}
			});
		}
	});
}
