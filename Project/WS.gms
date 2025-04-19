$eolCom #
set j indexy promenych /1,2,3/;
set i indexy obmedzeni /1,2,3,4/;
set s pocet scenaru /1,2/;  # TODO: Zatial nepouzitelne
parameter c(j) /1 90,2 150,3 155/; #Hodnota jednotkovych zisku 1 sivy riadok
parameter b(i) /1 9600,2 9600,3 10400,4 5400/; #Hodnota jednotlivych obmedzeni, purpurova hodnota
parameter p(s) /1 0.5, 2 0.5/; #Pravdepodobnost nastatnia nejakeho pripadu
parameter cs(j,s);
parameter as(i,j,s);
Parameter bs(i,s);
scalar EzWS, varzWS, szWS;
parameter zWSmax(s), xWSmax(j,s);
table a(i,j) #Dalsie riadky obmedzenia (pocet riadkov x pocet stlpcov)
    1 2  3
1  35 35 35    
2  14 13 15
3  38 35 34
4  25 20 20;
as(i,j,s) = uniform(0.8,2.0);
display "Test:";
display as;
Variable z;
Positive Variable x(j);
Equation ucelfce, omez(i);
ucelfce.. z =e= sum(j, c(j) * x(j));
omez(i).. sum(j, a(i,j)* x(j)) =l= b(i);
Model vyroba /ucelfce, omez/;
Solve vyroba maximizing z using LP;
display z.L, x.L;
File out /"Vysledky.txt"/;
put out;
put "Vysledky a vstupy:"/;#/ newline -> what the fuck is this language
put "------------------"/;
put "S ","p(s) ","opt? ","num? ", "z_max "/;
loop (j, put "x(",j.TL:1,") ";); #
loop (j, put "c(",j.TL:1,") ";);
loop (i, put "b(",i.TL:1,") ";);
loop (i, put "a(",i.TL:1,") ";);
loop (i, loop(j, put "a(", i.TL:1, ",", j.TL:1,") ";););
put /;