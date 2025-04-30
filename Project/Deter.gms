$eolCom #
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 4000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    c(j) /1 145, 2 150, 3 147/;
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
Variable z;
Positive Variables x(j);
Equations ucelfce, omez0(i), omez1(i), omez2(j), omez3(j);
ucelfce..           z =E= sum( j, c(j) * x(j) );
omez0(i)..          sum( j, a(i,j) * x(j) ) =L= b_UP(i);
omez1(i)..          sum( j, a(i,j) * x(j) ) =G= b_DOWN(i);
omez2(j)..          x(j) =L= x_UP(j);
omez3(j)..          x(j) =G= x_DOWN(j);
model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;
solve vyroba maximizing z using LP;
display z.L, x.L;
file out / "vysledky.txt" /;
put out;
put "Vysledky a vstupy:" /;
put "==================" / /;
put "opt?":7, "num?":7, "z_max":12;
loop(j, put "x(" , j.tl:0 , ")":10;);
loop(i, put "b_UP(" , i.tl:0 , ")":10;);
loop(i, put "b_DN(" , i.tl:0 , ")":10;);
loop(j, put "x_UP(" , j.tl:0 , ")":10;);
loop(j, put "x_DN(" , j.tl:0 , ")":10;);
loop(i, loop(j, put "a(" , i.tl:0 , "," , j.tl:0 , ")":10;););
put /;


put vyroba.modelstat:2:0, vyroba.solvestat:7:0, z.L:12:2
loop(j,
    put x.L(j):10:2;
    put "":1;
    put$(ord(j)=1) "":2;
    put$(ord(j)=2) "":2;
);
loop(i,
    put b_UP(i):14:2;
    put$(ord(i)=1) "":2;
    put$(ord(i)=2) "":2;
    put$(ord(i)=3) "":2;
);
loop(i, put b_DOWN(i):16:2;);
loop(j, put x_UP(j):16:2;);
loop(j, put x_DOWN(j):16:2;);
loop(i, loop(j, put a(i,j):15:2;););
put /;