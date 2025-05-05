$eolCom #
$onundf
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
set s / 1 * 3/;
scalar EzEOAS, varzEOAS, szEOAS;
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 4000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    qp(i) /1 0, 2 0, 3 0/,
    qm(i) /1 10, 2 20, 3 30/, #Naklady na zdroje navic -> TODO: Opytat sa teamu
    qms(i,s) "naklad na zdroje navic"  /1.1 10, 2.1 20, 3.1 30,
                                        1.2 20, 2.2 30, 3.2 10,
                                        1.3 40, 2.3 20, 3.3 10/,
    qps(i,s) "naklad na zdroje navic"  /1.1 0, 2.1 0, 3.1 0,
                                        1.2 0, 2.2 0, 3.2 0,
                                        1.3 0, 2.3 0, 3.3 0/,
    c(j) /1 145, 2 150, 3 147/
    p(s) /1 0.35,2 0.35,3 0.3/,
    as(i,j,s)
    /1.1.1 30, 1.2.1 28, 1.3.1 27, 
     2.1.1 22, 2.2.1 23, 2.3.1 21,
     3.1.1 31, 3.2.1 34, 3.3.1 34,
     4.1.1 24, 4.2.1 20, 4.3.1 24,
     1.1.2 30, 1.2.2 27, 1.3.2 26,
     2.1.2 22, 2.2.2 23, 2.3.2 21,
     3.1.2 31, 3.2.2 35, 3.3.2 34,
     4.1.2 23, 4.2.2 21, 4.3.2 22,
     1.1.3 25, 1.2.3 30, 1.3.3 30,
     2.1.3 25, 2.2.3 20, 2.3.3 22,
     3.1.3 30, 3.2.3 35, 3.3.3 32,
     4.1.3 18, 4.2.3 20,4.3.3 22/,
    zEOASmax(s),
    xEOASmax(j,s);
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
Variable z;
Positive Variables x(j),yp(i),ym(i),yps(i,s),yms(i,s);
x.up(j) = x_UP(j);
x.lo(j) = x_DOWN(j);
Equations ucelfce, ucelfceEO,omez1(i) ,omez2(i), omez1AS(i,s),omez2AS(i,s);

ucelfceEO..    z =E= sum(j,c(j)*x(j)) - sum(s,p(s)*sum(i,qms(i,s)*yms(i,s) + qps(i,s)*yps(i,s)));          
ucelfce..      z =E= sum(j,c(j)*x(j)) - sum(i,qm(i)*ym(i)+qp(i)*yp(i));
                          
omez1(i)..          sum(j, a(i,j) * x(j)) + yp(i) - ym(i) =E= b_UP(i);
omez2(i)..          sum(j, a(i,j) * x(j)) =G= b_DOWN(i);
omez1AS(i,s)..      sum(j, as(i,j,s) * x(j)) + yps(i,s) - yms(i,s) =E= b_UP(i);
omez2AS(i,s)..      sum(j, as(i,j,s) * x(j)) =G= b_DOWN(i);

model vyroba / ucelfceEO, omez1AS, omez2AS /;
model verifikace / ucelfce, omez1,omez2 /;
yp.UP(i) = 1000;
ym.UP(i) = 1000;
yps.UP(i,s) = 1000;
yms.UP(i,s) = 1000;

file out / "vysledkyEOAS_TS.html" /;
put out;
put "<head>";
put '<link rel="stylesheet" href="styles.css">';
put "</head>";
put "Vysledky a vstupy:<br>" /;
put "==================<br>" / /;
put "<table>";
put "<tr>";
put "<th>s</th>",
    "<th>p(s)</th>",
    "<th>opt?</th>",
    "<th>num?</th>",
    "<th>z_max</th>";
loop(j, put "<th>x(",j.tl:0,")</th>";);
loop(i, put "<th>b_UP(",i.tl:0,")</th>";);
loop(i, put "<th>b_DN(",i.tl:0,")</th>";);
loop(i, put "<th>Yp(",i.tl:0,")</th>";);
loop(i, put "<th>Ym(",i.tl:0,")</th>";);
loop(j, put "<th>x_UP(",j.tl:0,")</th>";);
put "</tr>";

solve vyroba maximizing z using LP;
display z.L, x.L;

If((vyroba.modelstat EQ 4) OR (vyroba.modelstat EQ 19), zEOASmax(s) = -INF;);
xEOASmax(j,s) = x.L(j);

put "<tr>";
put "<td>EOAS</td>",
    "<td></td>",
    "<td>"vyroba.modelstat"</td>",
    "<td>"vyroba.solvestat"</td>",
    "<td>"z.L"</td>";
loop(j, put "<td>"x.L(j)"</td>";);

loop(s,
    If((ord(s)=1),
        loop(i,put "<td>"b_UP(i)"</td>";);
        loop(i,put "<td>"b_DOWN(i)"</td>";);
        loop(i,put "<td>"yp.L(i)"</td>";);
        loop(i,put "<td>"ym.L(i)"</td>";);
        loop(j,put "<td>"x_UP(j)"</td>";);
        put "</tr>";
    else
        put "<tr>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        put "<td></td>";
        loop(i,put "<td>"b_UP(i)"</td>";);
        loop(i,put "<td>"b_DOWN(i)"</td>";);
        loop(i,put "<td>"yp.L(i)"</td>";);
        loop(i,put "<td>"ym.L(i)"</td>";);
        loop(j,put "<td>"x_UP(j)"</td>";);
        put "</tr>";
    );
);
put "</table>";
put "<br>";
put "<table>";
put "<tr>";
loop(j, put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.TL:1,")</th>";);
loop(i, loop(j, put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
loop(i, put "<th>qp(",i.TL:1,")</th>";);
loop(i, put "<th>qm(",i.TL:1,")</th>";);
put "</tr>";

loop(s,
        put "<tr>";
        loop(j,put "<td>"x_DOWN(j)"</td>";);
        loop(j,put "<td>"c(j)"</td>";);
        loop(i, loop(j, put "<td>"as(i,j,s)"</td>";););
        loop(i, put "<td>"qps(i,s)"</td>";);
        loop(i, put "<td>"qms(i,s)"</td>";);
        put "</tr>";
);    
put "</table>";
put "<br>";

x.Lo(j) = x.L(j);
x.Up(j) = x.L(j);

put "<table>";
put "<tr>";
put "<th>s</th>",
    "<th>p(s)</th>",
    "<th>opt?</th>",
    "<th>num?</th>",
    "<th>z_max</th>";
loop(j, put "<th>x(",j.tl:0,")</th>";);
loop(i, put "<th>b_UP(",i.tl:0,")</th>";);
loop(i, put "<th>b_DN(",i.tl:0,")</th>";);
loop(i, put "<th>Yp(",i.tl:0,")</th>";);
loop(i, put "<th>Ym(",i.tl:0,")</th>";);
loop(j, put "<th>x_UP(",j.tl:0,")</th>";);
put "</tr>";

loop(s,
    a(i,j) = as(i,j,s);
    qp(i) = qps(i,s);
    qm(i) = qms(i,s);
    
    solve verifikace maximizing z using LP;
    display z.L, x.L;
    zEOASmax(s) = z.L;
    
    If(vyroba.modelstat EQ 4, zEOASmax(s) = -INF;);
    xEOASmax(j,s) = x.L(j);
    put "<tr>"
    put "<td>"s.tl"</td>",
        "<td>"p(s)"</td>",
        "<td>"vyroba.modelstat"</td>",
        "<td>"vyroba.solvestat"</td>",
        "<td>"z.L"</td>";
    loop(j,put "<td>"x.L(j)"</td>";);
    loop(i,put "<td>"b_UP(i)"</td>";);
    loop(i,put "<td>"b_DOWN(i)"</td>";);
    loop(i,put "<td>"yp.L(i)"</td>";);
    loop(i,put "<td>"ym.L(i)"</td>";);
    loop(j,put "<td>"x_UP(j)"</td>";);
    put "</tr>";  
);
put "</table>";
put "<br>";


put "<table>";
put "<tr>";
loop(j, put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.TL:1,")</th>";);
loop(i, loop(j, put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
loop(i, put "<th>qp(",i.TL:1,")</th>";);
loop(i, put "<th>qm(",i.TL:1,")</th>";);
put "</tr>";

loop(s,
    a(i,j) = as(i,j,s);
    qp(i) = qps(i,s);
    qm(i) = qms(i,s);
    
    solve verifikace maximizing z using LP;
    display z.L, x.L;
    zEOASmax(s) = z.L;
    
    If(vyroba.modelstat EQ 4, zEOASmax(s) = -INF;);
    xEOASmax(j,s) = x.L(j);
    put "<tr>"
    loop(j,put "<td>"x_DOWN(j)"</td>";);
    loop(j,put "<td>"c(j)"</td>";);
    loop(i,loop(j, put "<td>"a(i,j)"</td>";););
    loop(i,put "<td>"qp(i)"</td>";);
    loop(i,put "<td>"qm(i)"</td>";);
    put "</tr>";  
);
put "</table>";
put "<br>";

EzEOAS   = sum(s,p(s) * zEOASmax(s) );
If(EzEOAS GT -INF,varzEOAS = sum(s,p(s) * zEOASmax(s) * zEOASmax(s) ) - EzEOAS * EzEOAS;);
If(EzEOAS EQ -INF,varzEOAS = UNDF;);
szEOAS   = sqrt(varzEOAS);

put "<table>";
put "<tr>";
put "<th>zTSmax</th>";
put "<th>varzEOAS</th>";
put "<th>sthzEOAS</th>";
put "</tr>";
put "<tr>";
put "<td>"EzEOAS"</td>";
put "<td>"varzEOAS"</td>";
put "<td>"szEOAS"</td>";
put "</tr>";
put "</table>";