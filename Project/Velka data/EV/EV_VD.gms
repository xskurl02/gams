$eolCom #
$onundf
set i /1 * 4/; #riadky x
set j /1 * 4/; #stlpce y
set s /1 * 3/; #scenare s
alias(s,s1);
scalar EzEV, varzEV, szEV;
parameter
    b_UP /1 4000,2 3000,3 4500,4 3000/, #Horna hranica suctu v riadku
    b_DOWN /1 1000,2 1500,3 3000,4 800/, #Dolna hranica suctu v riadku
    x_UP /1 100,2 100,3 100/, # Horna hranica x
    x_DOWN /1 20,2 40,3 20/, # Dolna hranica x
    p(s) /1 0.35,2 0.35,3 0.3/,
    c(j) /1 145, 2 150, 3 147, 4 143/,
    as(i,j,s),
    zEVmax(s),
    xEVmax(j,s);
table a(i,j)
    1   2   3   4
1   30  28  27  20
2   22  23  21  21
3   31  34  34  23
4   24  20  24  20;
*-------------------------
*Large Data Creation
*-------------------------
as(i,j,s) = uniform(15,30);
p(s) = uniform(0,1);
p(s) = p(s)/sum(s1,p(s1));
*--------------------------
display as;
Variable z;
Positive Variables x(j);
Equations ucelfce, omez0(i), omez1(i), omez2(j), omez3(j);
ucelfce..           z =E= sum(j,c(j) * x(j));
omez0(i)..          sum(j,a(i,j) * x(j)) =L= b_UP(i);
omez1(i)..          sum(j,a(i,j) * x(j)) =G= b_DOWN(i);
omez2(j)..          x(j) =L= x_UP(j);
omez3(j)..          x(j) =G= x_DOWN(j);
model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;

file out / "vysledky.txt" /;
put out;
put "Vysledky a vstupy:" /;
put "==================" / /;
put "s":3, "p(s)":6, "opt?":7, "num?":7, "z_max":12;
loop(j, put "x(" , j.tl:0 , ")":10;);
loop(i, put "b_UP(" , i.tl:0 , ")":10;);
loop(i, put "b_DN(" , i.tl:0 , ")":10;);
loop(j, put "x_UP(" , j.tl:0 , ")":10;);
loop(j, put "x_DN(" , j.tl:0 , ")":10;);
loop(i, loop(j, put "a(" , i.tl:0 , "," , j.tl:0 , ")":10;););
put /;

a(i,j) = sum(s, p(s) * as(i,j,s));

solve vyroba maximizing z using LP;
display z.L, x.L;

* ---------- Expected-value row (base LP) ----------
put "EV":3, "":1, vyroba.modelstat:7:0, vyroba.solvestat:7:0, z.l:12:2;
loop(j,
    put x.l(j):10:2;
    put "":2
    put$(ord(j)=1) "":1;
    put$(ord(j)=2) "":1;
);

loop(i,
    put b_UP(i):13:2;
    put$(ord(i)=1) "":3;
    put$(ord(i)=2) "":3;
    put$(ord(i)=3) "":3;
);

loop(i,put b_DOWN(i):16:2;);
loop(j,put x_UP(j):16:2;);
loop(j,put x_DOWN(j):16:2;);
loop(i,loop(j,put a(i,j):15:2;););
put /;

#To Find the current value of X and cap the fucker
x.Lo(j) = x.L(j);
x.Up(j) = x.L(j);

loop(s,
    a(i,j) = as(i,j,s)
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    zEVmax(s) = z.L;
    If((vyroba.modelstat EQ 4) OR (vyroba.modelstat EQ 19), zEVmax(s) = -INF;);
    put s.tl:3, p(s):5:3, vyroba.modelstat:3:0, vyroba.solvestat:7:0, z.L:11:2;
    loop(j,
        put x.L(j):11:2;
        put "":1;
        put$(ord(j)=1) "":1;
        put$(ord(j)=2) "":1;
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
);

* ---------- risk measures ----------
EzEV   = sum(s,p(s) * zEVmax(s) );
put "EzEV   = ", EzEV:12:2 /;
varzEV = sum(s,p(s) * zEVmax(s) * zEVmax(s) ) - EzEV * EzEV;
put "varzEV = ", varzEV:12:2 /;
szEV   = sqrt(varzEV);
put "stdzEV = ", szEV:12:2 /;