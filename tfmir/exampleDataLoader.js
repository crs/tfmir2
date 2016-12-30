function wait(ms){
   var start = new Date().getTime();
   var end = start;
   console.log("Waiting");
   while(end < start + ms) {
     end = new Date().getTime();
  }
   console.log("Continue");
}

function loadExampleData() {
	console.log('Mimic sample data upload');
	$('#mirnaDemo').attr('value', true);
	$('#mrnaDemo').attr('value', true);

	uploadFile('miRNA');
	wait(1500);
	uploadFile('mRNA');	

	$('#mirnaDemo').attr('value', false);
	$('#mrnaDemo').attr('value', false);

	$('select[id=disease] option[selected="selected"]').removeAttr('selected');
	$('select[id=tissue] option[selected="selected"]').removeAttr('selected');
	$('select[id=function] option[selected="selected"]').removeAttr('selected');
	//$('select[id=disease] option[value="Breast Neoplasms"]').attr('selected','selected');
	
	var text1 = 'Breast Neoplasms';
	var text2 = 'Breast - Mammary Tissue';
	var text3 = 'cell proliferation';


	$("select[id=disease] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');

	$("select[id=disease] option").filter(function() {
    //may want to use $.trim in here
    	return $(this).text() == text1; 
	}).attr('selected', true).addClass('highlight');

	//$('#disease').chosen().change();
	$("#disease").trigger("chosen:updated");
	$("select[id=disease] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');
	
	$("select[id=tissue] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');
	
	$("select[id=tissue] option").filter(function() {
	//may want to use $.trim in here
		 return $(this).text() == text2; 
	}).attr('selected', true).addClass('highlight');

	$("#tissue").trigger("chosen:updated");
	$("select[id=tissue] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');
	
	$("select[id=function] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');

	$("select[id=function] option").filter(function() {
		 //may want to use $.trim in here
		 return $(this).text() == text3; 
	 }).attr('selected', true).addClass('highlight');

	//$('#disease').chosen().change();
	$("#function").trigger("chosen:updated");
	$("select[id=function] option").filter(function()
	{
		return $(this).hasClass('highlight');
	}).removeClass('highlight').removeAttr('selected');

	$("select[id=evidence] option").filter(function() {
    //may want to use $.trim in here
    	return $(this).text() == 'Experimental'; 
	}).attr('selected', true).addClass('highlight');
	
	$("select[id=species] option").filter(function() {
	 //may want to use $.trim in here
	    return $(this).text() == 'Human'; 
	}).attr('selected', true).addClass('highlight');
	//$('#disease').chosen().change();
	
	$("#evidence").trigger("chosen:updated");
	$("#species").trigger("chosen:updated");
	$('#disease').val("Breast Neoplasms");
	$('#tissue').val("Breast - Mammary Tissue");
	$('#function').val("cell proliferation");
	$('#evidence').val("Experimental");
	$('#species').val("Human");
	$('#enrich_pvalue').val("0.05");
	$('#ppIcut').val("0.8");
	$('#orapvalue').val("0.05");
	


	docCookies.setItem('evidence', 'Experimental');
	docCookies.setItem('species', 'Human');
	docCookies.setItem('randomization.method', 'conserved');
	docCookies.setItem('disease', 'Breast Neoplasms');
	docCookies.setItem('tissue', 'Breast - Mammary Tissue');
	docCookies.setItem('function', 'cell proliferation');
	docCookies.setItem('enrich_pvalue', '0.05');
	docCookies.setItem('ppIcut', '0.8');
	docCookies.setItem('orapvalue', '0.05');
	log("Example files for breast cancer have been loaded into your session. Check your p-value and experimental evidence and click the processing button to start analysis.")
}
