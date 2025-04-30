$eolCom #
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 4000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    qm(i) / 1 10, 2 20, 3 30/, #Naklady na zdroje navic -> TODO: Opytat sa teamu
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
ucelfce..           z =E= sum(j, c(j) * x(j)) - sum(i, qm(i)* ym(i) );
omez0(i)..          sum(j,a(i,j) * x(j)) + yp(i) - ym(i) =E= b_UP(i);
omez1(i)..          sum(j,a(i,j) * x(j)) =G= b_DOWN(i);
omez2(j)..          x(j) =L= x_UP(j);
omez3(j)..          x(j) =G= x_DOWN(j);
model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;
yp.UP(i) = 100;
ym.UP(i) = 100;
solve vyroba maximizing z using LP;
display z.L, x.L;
file out / "vysledky.txt" /;
put out;
put "Vysledky a vstupy:" /;
put "==================" / /;
put 'opt?':7,'num?':7,'z_max':12;
loop(j, put "x(", j.TL:1, ")":10;);
loop(i, put "yp(", i.TL:1, ")":10;);
loop(i, put "ym(", i.TL:1, ")":10;);
loop(j, put "c(", j.TL:1, ")":10;);
loop(i, put "b_up(", i.TL:1, ")":10;);
loop(i, put "b_down(", i.TL:1, ")":10;);
loop(i,loop(j, put " a(",i.TL:1,",",j.TL:1, ")":10;););
put /;
put vyroba.modelstat:2:0,vyroba.solvestat:7:0,z.l:12:2;

loop(j,
    put x.L(j):10:2;
    put "":1;
    put$(ord(j)=1) "":2;
    put$(ord(j)=2) "":2;
);
loop(i,
    put yp.l(i):11:2;
    put "":1;
    put$(ord(i)=1) "":2;
    put$(ord(i)=2) "":2;
    put$(ord(i)=3) "":2;
);
loop(i, put ym.l(i):14:2;);
loop(j, put c(j):13:2;);
loop(i,
    put b_up(i):13:0;
    put$(ord(i)=1) "":3;
    put$(ord(i)=2) "":3;
    put$(ord(i)=3) "":3;
);
loop(i,
    put b_down(i):16:0 ;
    put$(ord(i)=1) "":2;
    put$(ord(i)=2) "":2;
    put$(ord(i)=3) "":1;
);
loop(i,
     loop(j,
        put$(ord(i)=1) "":1;
        put a(i,j):16:0;
        );   
);
