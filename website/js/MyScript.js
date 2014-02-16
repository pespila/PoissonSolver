if (navigator.appName == 'Microsoft Internet Explorer') {
	document.createElement("nav");
	document.createElement("article");
	document.createElement("section");
	document.createElement("aside");
	document.createElement("header");
	document.createElement("footer");
}

var counter=0;
$(document).ready(function(){
    $("#maincontent").html("This is an online solver for the disretized poisson equation on the unity square (2D). Click \"Run\" to solve an equation. On the right side you'll find some links which display the used algorithms. To download source code or run this programme under Linux visit my <a href=\"http://www.github.com/BauerMichael/PoissonSolver\" target=\"_blank\">Github Repository</a>");
});

$("#home").click(function(){
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").text("PoissonSolver");
    $("#maincontent").html("This is an online solver for the disretized poisson equation on the unity square (2D). Click \"Run\" to solve an equation. On the right side you'll find some links which display the used algorithms. To download source code or run this programme under Linux visit my <a href=\"http://www.github.com/BauerMichael/PoissonSolver\" target=\"_blank\">Github Repository</a>");
});
$("#about").click(function(){
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").text("Bachelor Thesis");
    $("#maincontent").html("<b>Michael Bauer</b>, degree course: <b>Computational Science</b>.<br><br>This PoissonSolver which has come to my mind while writing and programming on my Bachelor Thesis with the title <b><i>Iterative Methods for solving the discretized poisson equation</b></i>. This thesis was written under the supervisor Prof. Dr. Harald Garcke.");
});
$("#contact").click(function(){
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").text("Contact me");
    $("#maincontent").html("E-Mail: <a href=mailto:michael1.bauer@stud.uni-r.de>Michael Bauer</a>");
});

$("#equ").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Possion Equation");
    $("#maincontent").html( '<img src="./img/poissonEqu.png" width=100% />' );
});
$("#jaco").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Jacobi Method");
    $("#maincontent").html( '<img src="./img/jacobi.png" width=100% />' );
});
$("#cg").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Conjugate Gradient");
    $("#maincontent").html( '<img src="./img/cg.png" width=100% />' );
});
$("#pcg").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Preconditioned Conjugate Gradient");
    $("#maincontent").html( '<img src="./img/pcg.png" width=100% />' );
});
$("#icg").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Incomplete Cholesky Decomposition");
    $("#maincontent").html( '<img src="./img/icg.png" width=100% />' );
});
$("#micg").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("Modified Incomplete Cholesky Decomposition");
    $("#maincontent").html( '<img src="./img/micg.png" width=100% />' );
});
$("#vc").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("V-cycle");
    $("#maincontent").html( '<img src="./img/vcycle.png" width=60% />' );
});
$("#wc").click(function() {
    $(".form").css({"visibility":"hidden","height":"0px"});
    $("#run").css({"visibility":"visible","height":"auto"});
    $("#calc").css({"visibility":"hidden","height":"0px","width":"0px"});
    $("h1").html("W-cycle");
    $("#maincontent").html( '<img src="./img/wcycle.png" width=55% />' );
});

$("#run").click(function() {
    $("h1").html("PoissonSolver");
    $("#maincontent").html("");
    $(".form").css({"visibility":"visible","height":"auto"});
    $("#run").css({"visibility":"hidden","height":"0px"});
    $("#calc").css({"visibility":"visible","height":"auto","width":"auto"});
});

$(document).ready(function () {

    function exportTableToCSV($table, filename) {

        var $rows = $table.find('tr:has(td)'),

            // Temporary delimiter characters unlikely to be typed by keyboard
            // This is to avoid accidentally splitting the actual contents
            tmpColDelim = String.fromCharCode(11), // vertical tab character
            tmpRowDelim = String.fromCharCode(0), // null character

            // actual delimiter characters for CSV format
            colDelim = '","',
            rowDelim = '"\r\n"',

            // Grab text from table into CSV formatted string
            csv = '"' + $rows.map(function (i, row) {
                var $row = $(row),
                    $cols = $row.find('td');

                return $cols.map(function (j, col) {
                    var $col = $(col),
                        text = $col.text();

                    return text.replace('"', '""'); // escape double quotes

                }).get().join(tmpColDelim);

            }).get().join(tmpRowDelim)
                .split(tmpRowDelim).join(rowDelim)
                .split(tmpColDelim).join(colDelim) + '"',

            // Data URI
            csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

        $(this)
            .attr({
            'download': filename,
                'href': csvData,
                'target': '_blank'
        });
    }

    // This must be a hyperlink
    $("#clickExcel").on('click', function (event) {
        alert("After opening, be sure that your selcted character set ist UTF-8!");
        // CSV
        exportTableToCSV.apply(this, [$('#myTable>table'), 'export.csv']);
        
        // IF CSV, don't do event.preventDefault() or return false
        // We actually need this to be a typical hyperlink
    });
});

$(document).ready(function() {
    $("#grid").keyup(function() {
        var x=$("#grid").val();
        var gridWidth=["2","4","8","16","32","64","128","256","512","1024","2048"];
        var proof=0;
        for(var i=0;i<gridWidth.length;i++) {
            if(x==gridWidth[i]) proof=1;
        }
        if(proof==0) {
            $("#opt [value=4],[value=5],[value=6]").attr("disabled","disabled");
        }
        if(proof==1) {
            $("#opt option").each(function() {
                $(this).attr("disabled",false);
            });
        }
        if(x>128) {
            $("#opt [value=7],[value=8],[value=9]").attr("disabled","disabled");
        }
        if(x>256) {
            $("#opt [value=10]").attr("disabled","disabled");
        }
        if(x>320) {
            $("#opt [value=1]").attr("disabled","disabled");
        }
    });
});

$(document).ready(function() {
    $("#calc").click(function() {
        $("#wait").css("visibility","visible");
        var x=$("#grid").val();
        var x=x-1;
        var y=$("#opt").val();
        var z=$(":input:checked").val();
        if(z !== undefined && y != 99) {
            $.ajax({
                type: "GET",
                url: "./cgi/poisson.cgi",
                data: { arg: x, alg: y, func: z }
            })
            .done(function( data,status ) {
                if(status == "success") {
                    counter++;
                    for(var i=1;i<counter;i++) {
                        $("#"+i).html("<button class=\"btn btn-danger btn-sm\">overwritten")
                    }
                    $(".table").append("<tr><td>"+counter+"</td>"+data+"<td id="+counter+"><a href=\"./Plot/plot.dat\" class=\"btn btn-primary btn-sm\" target=_blank>data file</a></td></tr>");
                    $("#wait").css("visibility","hidden");
                } else {
                    alert("Couldn't run programme!");
                }
            });
        } else if(z !== undefined && y == 99) {
            var myArray=[];
            var u=0;
            var wait=0;
            if(x<=128) {
                element=$("#opt :enabled");
                element.each(function() {
                    if($(this).val()!=99) {
                        myArray[u]=$(this).val();
                        u++;
                    }
                });
            }
            for(var i=0;i<myArray.length;i++) {
                $.ajax({
                    type: "GET",
                    url: "./cgi/poisson.cgi",
                    data: { arg: x, alg: myArray[i], func: z }
                })
                .done(function( data,status ) {
                    if(status == "success") {
                        counter++;
                        for(var i=1;i<counter;i++) {
                            $("#"+i).html("<button class=\"btn btn-danger btn-sm\">overwritten")
                        }
                        $(".table").append("<tr><td>"+counter+"</td>"+data+"<td id="+counter+"><a href=\"./Plot/plot.dat\" class=\"btn btn-primary btn-sm\" target=_blank>data file</a></td></tr>");
                        wait++;
                        if(wait==myArray.length) $("#wait").css("visibility","hidden");
                    } else {
                        alert("Couldn't run programme!");
                    }
                });
            }
        } else {
            alert("Sorry, missing input!");
        }
    });
});