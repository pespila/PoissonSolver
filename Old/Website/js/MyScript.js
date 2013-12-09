if (navigator.appName == 'Microsoft Internet Explorer') {
	document.createElement("nav");
	document.createElement("article");
	document.createElement("section");
	document.createElement("aside");
	document.createElement("header");
	document.createElement("footer");
}

$(document).ready(function(){
	$("#maincontent").text("With this website you can solve the Discretized Poisson Equation with the help of the Conjugate Gradient Method.");
});
$("#home").click(function(){
	$("h1").text("PoissonSolver Online");
	$("#maincontent").text("With this website you can solve the Discretized Poisson Equation with the help of the Conjugate Gradient Method.");
});
$("#about").click(function(){
	$("h1").text("Bachelor Thesis of...");
	$("#maincontent").text("...Michael Bauer, degree course: Computational Science, title: ...");
});
$("#contact").click(function(){
	$("h1").text("Contact me");
	$("#maincontent").html("E-Mail: <a href=mailto:michael1.bauer@stud.uni-r.de>Michael Bauer</a>");
});
$("#run").click(function(){
	
});