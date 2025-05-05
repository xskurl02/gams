$eolCom #
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 4000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    qm(i) / 1 10, 2 20, 3 30, 4 40/, #Naklady na zdroje navic -> TODO: Opytat sa teamu
    c(j) /1 145, 2 150, 3 147/;
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
Variable z;
Positive Variables x(j),yp(i),ym(i);
Equations ucelfce, omez0(i), omez1(i), omez2(j), omez3(j);
ucelfce..   z =E= sum(j, c(j) * x(j)) - sum(i, qm(i)* ym(i));
omez0(i)..  sum(j,a(i,j) * x(j)) + yp(i) - ym(i) =E= b_UP(i);
omez1(i)..  sum(j,a(i,j) * x(j)) =G= b_DOWN(i);
omez2(j)..  x(j) =L= x_UP(j);
omez3(j)..  x(j) =G= x_DOWN(j);

model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;
yp.UP(i) = 100;
ym.UP(i) = 100;

solve vyroba maximizing z using LP;
display z.L, x.L;

file out / "vysledkyDeter_TS.html" /;
put out;
put "<head>";
put '<link rel="stylesheet" href="styles.css">';
put "</head>";

put "Vysledky a vstupy:<br>" /;
put "==================<br>" /;
put "<table>";
put "<tr>";
put '<th>opt?</th>',
    '<th>num?</th>',
    '<th>z_max</th>';
loop(j,put "<th>x(",j.tl:0,")</th>";);
loop(i,put "<th>yp(",i.tl:0,")</th>";);
loop(i,put "<th>ym(",i.tl:0,")</th>";);
loop(i,put "<th>b_UP(",i.tl:0,")</th>";);
loop(i,put "<th>b_DN(",i.tl:0,")</th>";);
put "</tr>";

put "<tr>";
put "<td>"vyroba.modelstat"</td>",
    "<td>"vyroba.solvestat"</td>",
    "<td>"z.l"</td>";
loop(j,put "<td>"x.L(j)"</td>";);
loop(i,put "<td>"yp.l(i)"</td>";);
loop(i,put "<td>"ym.l(i)"</td>";);
loop(i,put "<td>"b_UP(i)"</td>";);
loop(i,put "<td>"b_DOWN(i)"</td>";);
put "</tr>";
put "</table>";
put "<br>";

put "<table>";
put "<tr>";
loop(j,put "<th>x_UP(",j.tl:0,")</th>";);
loop(j,put "<th>x_DN(",j.tl:0,")</th>";);
loop(j,put "<th>c(",j.TL:1,")</th>";);
loop(i,loop(j, put "<th>a(",i.tl:0,",",j.tl:0,")</th>";););
put "</tr>";

put "<tr>";
loop(j,put "<td>"x_UP(j)"</td>";);
loop(j,put "<td>"x_DOWN(j)"</td>";);
loop(j,put "<td>"c(j)"</td>";);
loop(i,loop(j,put "<td>"a(i,j)"</td>";););
put "</tr>";
put "</table>";

put "<br>";
put "<table>";
put "<tr>";
loop(i, put "<th>qm(",i.TL:1,")</th>";);
put "</tr>";
put "<tr>";
loop(i,put "<td>"qm(i)"</td>";);
put "</tr>";
put "</table>";