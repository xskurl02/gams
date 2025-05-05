$eolCom #
set i /1 * 4/; #stlpce x
set j /1 * 4/; #riadky y
set s /1 * 3/; #scenare s
alias(s,s1);
scalar EzWS, varzWS, szWS;
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 2500,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100,4 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20,4 20/, # Dolna hranica x
    p(s) /1 0.35,2 0.35,3 0.3/, # Pravdepodobnost scernaru
    c(j) /1 145,2 150,3 147,4 143/, #Cena Kavy
    as(i,j,s),
    zWSmax(s),
    xWSmax(j,s);
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
*-------------------------
*Large Data Creation
*-------------------------
as(i,j,s) = uniform(10,30);
p(s) = uniform(0,1);
p(s) = p(s)/sum(s1,p(s1));
*--------------------------
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


*--------------------------
* OUTPUT Formatting
*--------------------------
file out / "vysledky_WS_VD.html" /;
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
put "</tr>";

loop(s,
    a(i,j) = as(i,j,s);
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    zWSmax(s) = z.L;
    xWSmax(j,s) = x.L(j);
    put "<tr>";
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
    put "<tr>";
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