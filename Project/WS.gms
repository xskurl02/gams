$eolCom #
set j indexy promenych /1,2,3/;
set i indexy obmedzeni /1,2,3,4/;
set s pocet scenaru /1,2,3/;  # TODO: Zatial nepouzitelne
alias(s,s1);
parameter c(j) /1 90,2 150,3 155/; #Hodnota jednotkovych zisku 1 sivy riadok
parameter b(i) /1 9600,2 9600,3 5400,4 5400/; #Hodnota jednotlivych obmedzeni, purpurova hodnota
parameter p(s) /1 0.5, 2 0.2, 3 0.3/; #Pravdepodobnost nastatnia nejakeho pripadu
Parameter dm(j) dolni mez        /1 20, 2 40, 3 20/,
          hm(j) horni mez        /1 170, 2 170, 3 170/;
parameter cs(j,s)/1.1 40, 1.2 20,1.3 20, 2.1 20, 2.2 40, 2.3 20, 3.1 20, 3.2 20, 3.3 40/;
parameter as(i,j,s);
Parameter bs(i,s)/1.1 9600, 2.1 9600, 3.1 5400, 4.1 5400, 1.2 9600, 2.2 9600, 3.2 5400, 4.2 5400, 1.3 6900, 2.3 6900, 3.3 5400, 4.3 5400/;
scalar EzWS, varzWS, szWS;
parameter zWSmax(s), xWSmax(j,s);
table a(i,j) #Dalsie riadky obmedzenia (pocet riadkov x pocet stlpcov)
    1 2  3
1  35 35 35    
2  14 13 15
3  38 35 34
4  25 20 20;
as(i,j,s) = uniform(0.8,2.0);
cs(j,s) = uniform(20,40);
bs(i,s)= uniform(5400,9600);
p(s) = uniform(0,1);
p(s) = p(s)/sum(s1,p(s1));
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
put "s","  p(s)"," opt?"," num?", "  z_max";
loop (j, put "  x(",j.TL:1,")";); #
loop (j, put "  c(",j.TL:1,")";);
loop (i, put "  b(",i.TL:1,")";);
loop (i, loop(j, put " a(", i.TL:1,",", j.TL:1,")"  ;););
put /;
loop(
    s, c(j) = cs(j,s);
    b(i) = bs(i,s);
    a(i,j) = as(i,j,s);
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    zWSmax(s) = z.L;
    xWSmax(j,s) = x.L(j);
    put s.TL:1"  ", p(s):4:0"  ", vyroba.modelstat:5:0"  ", vyroba.solvestat:5:0"  ",
    z.L:7:1; loop(j, put x.L(j):6:0;);
    loop(j, put c(j):6:0;);
    loop(i, put b(i):7:0;); 
    loop(i, loop(j,put a(i,j):8:0;););
    put /;
);
put /; EzWS =sum(s,p(s) * zWSmax(s) );
put "EzWS   = ", EzWS:7:0 /;
varzWS = sum(s,p(s) * zWSmax(s) * zWSmax(s) ) - EzWS * EzWS;
put "varzWS = ", varzWS:7:0 /;
szWS   = sqrt(varzWS);
put "stdzWS = ", szWS:7:0 /;