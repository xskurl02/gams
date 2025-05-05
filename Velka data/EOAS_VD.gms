$eolCom #
$onundf
set i /1 * 4/; #stlpce x
set j /1 * 4/; #riadky y
set s /1 * 3/; #scenare s
alias(s,s1);
scalar EzEOAS, varzEOAS, szEOAS;
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 2500,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100,4 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20,4 20/, # Dolna hranica x
    p(s) /1 0.35,2 0.35,3 0.3/,
    c(j) /1 145, 2 150, 3 147, 4 143/,
    as(i,j,s),
    zEOASmax(s),
    xEOASmax(j,s);
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   26  20  24;

*-------------------------
*Large Data Creation
*-------------------------
as(i,j,s) = uniform(15,30);
p(s) = uniform(0,1);
p(s) = p(s)/sum(s1,p(s1));

Variable z;
Positive Variables x(j);
*-------------------------------- 
* upper / lower bounds for x
*-------------------------------- 
x.up(j) = x_UP(j);
x.lo(j) = x_DOWN(j);
Equations ucelfce, ucelfceEO,omez1(i) ,omez2(i), omez1AS(i,s),omez2AS(i,s);
ucelfceEO..   z =E=    sum(s, p(s) * sum( j,    c(j)   * x(j) ));
ucelfce..     z =E=                  sum( j,    c(j)   * x(j) );
omez1(i)..             sum(j,    a(i,j) * x(j) ) =L= b_UP(i);
omez2(i)..             sum(j,    a(i,j) * x(j) ) =G= b_DOWN(i);
omez1AS(i,s)..         sum(j, as(i,j,s) * x(j) ) =L= b_UP(i);
omez2AS(i,s)..         sum(j, as(i,j,s) * x(j) ) =G= b_DOWN(i);

*model vyroba / ucelfceEO, omez1AS, omez2AS /;
*model verifikace / ucelfce, omez1,omez2 /;

model vyroba / ucelfceEO, omez1AS/;
model verifikace / ucelfce, omez1/;
*---------------------------------------------------------------------------
*  >>> neatly formatted output file  <<<
*---------------------------------------------------------------------------
file out / "vysledky_EOAS_VD.html" /;
put out;
put "<head>";
put '<link rel="stylesheet" href="styles.css">';
put "</head>";
* ---------- column headings ----------
put "Vysledky a vstupy:<br>" /;
put "==================<br>" /;
put '<table id="EOAS1/2">';
put "<tr>";
put "<th>s</th>",
    "<th>p(s)</th>",
    "<th>opt?</th>",
    "<th>num?</th>",
    "<th>z_max</th>";
loop(j, put "<th>x(",j.tl:0,")</th>";);
loop(i, put "<th>b_UP(",i.tl:0,")</th>";);
loop(i, put "<th>b_DN(",i.tl:0,")</th>";);
loop(j, put "<th>x_UP(",j.tl:0,")</th>";);
put "</tr>"
#a(i,j) = sum(s, p(s) * as(i,j,s))
solve vyroba maximizing z using LP;
display z.L, x.L;
* ---------- Overall EOAS output(1/2) ----------
put "<tr>";
put "<td>EOAS</td>",
    "<td></td>",
    "<td>"vyroba.modelstat"</td>",
    "<td>"vyroba.solvestat"</td>",
    "<td>"z.L"</td>";
loop(j,put "<td>"x.L(j)"</td>";);
loop(i,put "<td>"b_UP(i)"</td>";);
loop(i,put "<td>"b_DOWN(i)"</td>";);
loop(j,put "<td>"x_UP(j)"</td>";);
put "</tr>";
put '</table id="EOAS1/2">';
put "<br>"
* ---------- column headings ----------
put '<table id="EOAS2/2">';
put "<tr>";
loop(j, put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.tl:0,")</th>";);
loop(i, loop(j, put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
put "</tr>";
* ---------- Overall EOAS output(2/2) ----------
put "<tr>";
loop(j,put "<td>"x_DOWN(j)"</td>";);
loop(j,put "<td>"c(j)"</td>";);
loop(i,loop(j,put "<td>"a(i,j)"</td>";););
put "</tr>";
put '</table id="EOAS2/2">';

#To Find the current value of X and cap the fucker
x.Lo(j) = x.L(j);
x.Up(j) = x.L(j);

* ---------- column headings ----------
put '<table id="Scenarios">';
put "<tr>";
put "<th>s</th>",
    "<th>p(s)</th>",
    "<th>opt?</th>",
    "<th>num?</th>",
    "<th>z_max</th>";
loop(j, put "<th>x(",j.tl:0,")</th>";);
loop(i, put "<th>b_UP(",i.tl:0,")</th>";);
loop(i, put "<th>b_DN(",i.tl:0,")</th>";);
loop(j, put "<th>x_UP(",j.tl:0,")</th>";);
put "</tr>"
put "<br>"
* ---------- Overall output(1/2) ----------
loop(s,
    a(i,j) = as(i,j,s)
    solve verifikace maximizing z using LP;
    display z.L, x.L;
    zEOASmax(s) = z.L;
    if ((vyroba.modelstat = 4) or (vyroba.modelstat = 19),zEOASmax(s) = -INF);
    put "<tr>"
    put "<td>"s.tl"</td>",
        "<td>"p(s)"</td>",
        "<td>"vyroba.modelstat"</td>",
        "<td>"vyroba.solvestat"</td>",
        "<td>"zEOASmax(s)"</td>";
    loop(j,put "<td>"x.l(j)"</td>";);
    loop(i,put "<td>"b_UP(i)"</td>";);
    loop(i,put "<td>"b_DOWN(i)"</td>";);
    loop(j,put "<td>"x_UP(j)"</td>";);
);
put '</table id="Scenarios">';
put "<br>";
* ------------ column headings -------------
put '<table id="Scenarios2">';
put "<tr>";
loop(j, put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.tl:0,")</th>";);
loop(i, loop(j, put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
put "</tr>";

loop(s,
    a(i,j) = as(i,j,s)
    solve verifikace maximizing z using LP;
    display z.L, x.L;
    zEOASmax(s) = z.L;
    if ((vyroba.modelstat = 4) or (vyroba.modelstat = 19),zEOASmax(s) = -INF);
    put "<tr>";
    loop(j,put "<td>"x_DOWN(j)"</td>";);
    loop(j,put "<td>"c(j)"</td>";);
    loop(i, loop(j, put "<td>"a(i,j)"</td>";););
    put "</tr>";
);
put '</table id="Scenarios2">';
put "<br>";

* ---------- risk measures ----------
EzEOAS=sum(s,p(s)*zEOASmax(s));
If(EzEOAS GT -INF,varzEOAS = sum(s,p(s) * zEOASmax(s) * zEOASmax(s) ) - EzEOAS * EzEOAS;);
If(EzEOAS EQ -INF,varzEOAS = UNDF;);
szEOAS   = sqrt(varzEOAS);

put "<table>";
put "<tr>";
put "<th>EzEOAS</th>";
put "<th>varzEOAS</th>";
put "<th>stdzEOAS</th>";
put "</tr>";
put "<tr>";
put "<td>"EzEOAS"</td>";
put "<td>"varzEOAS"</td>";
put "<td>"szEOAS"</td>";
put "</tr>";
put "</table>";