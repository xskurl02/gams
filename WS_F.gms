$eolCom #
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
set s /1 * 3/; #scenare s
scalar EzWS, varzWS, szWS;
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 4000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    p(s) /1 0.35,2 0.35,3 0.3/,
    c(j) /1 145, 2 150, 3 147/,
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
    4.1.3 18, 4.2.3 20, 4.3.3 22/,
    zWSmax(s),
    xWSmax(j,s);
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
Variable z;
Positive Variables x(j);
Equations ucelfce, omez0(i), omez1(i), omez2(j), omez3(j);
ucelfce..   z =E= sum( j, c(j) * x(j) );
omez0(i)..  sum( j, a(i,j) * x(j) ) =L= b_UP(i);
omez1(i)..  sum( j, a(i,j) * x(j) ) =G= b_DOWN(i);
omez2(j)..  x(j) =L= x_UP(j);
omez3(j)..  x(j) =G= x_DOWN(j);
model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;

solve vyroba maximizing z using LP;
display z.L, x.L;

file out / "vysledky_WS.html" /;
put out;
put "<head>";
put '<link rel="stylesheet" href="styles.css">';
put "</head>";

put "Vysledky a vstupy:<br>"/;
put "==================" / /;
put '<table id="Scenarios">';
put "<tr>";
put "<th>s</th>",
    "<th>p(s)</th>",
    "<th>opt?</th>",
    "<th>num?</th>",
    "<th>z_max</th>";
loop(j,put "<th>x(",j.tl:0,")</th>";);
loop(i,put "<th>b_UP(",i.tl:0,")</th>";);
loop(i,put "<th>b_DN(",i.tl:0,")</th>";);
loop(j,put "<th>x_UP(",j.tl:0,")</th>";);
put "</tr>"
put ""

loop(s,
    a(i,j) = as(i,j,s);
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    zWSmax(s) = z.L;
    xWSmax(j,s) = x.L(j);
    put "<tr>"
    put "<td>"s.tl"</td>",
        "<td>"p(s)"</td>",
        "<td>"vyroba.modelstat"</td>",
        "<td>"vyroba.solvestat"</td>",
        "<td>"z.L"</td>";
    loop(j,put "<td>"x.L(j)"</td>";);
    loop(i,put "<td>"b_UP(i)"</td>";);
    loop(i,put "<td>"b_DOWN(i)"</td>";);
    loop(j,put "<td>"x_UP(j)"</td>";);
    put "</tr>";
);
put "</table>";
put "<br>";

put '<table id="Scenarios2">';
put "<tr>";
loop(j,put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.tl:0,")</th>";);
loop(i,loop(j,put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
put "</tr>";
loop(s,
    a(i,j) = as(i,j,s);
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    put "<tr>"
    loop(j,put "<td>"x_DOWN(j)"</td>";);
    loop(j,put "<td>"c(j)"</td>";);
    loop(i,loop(j,put "<td>"a(i,j)"</td>";););
    put "</tr>";
);
put "</table>";
put "<br>";
EzWS=sum(s,p(s)*zWSmax(s));
varzWS=sum(s,p(s)*zWSmax(s)*zWSmax(s))-EzWS*EzWS;
szWS=sqrt(varzWS);

put "<table>";
put "<tr>";
put "<th>EzWS</th>";
put "<th>varzWS</th>";
put "<th>stdzWS</th>";
put "</tr>";
put "<tr>";
put "<td>"EzWS"</td>";
put "<td>"varzWS"</td>";
put "<td>"szWS"</td>";
put "</tr>";
put "</table>";