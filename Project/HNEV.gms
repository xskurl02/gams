$eolCom #
$onundf
set i /1 * 4/; #stlpce x
set j /1 * 3/; #riadky y
set s /1 * 3/; #scenare s
scalar EzEV, varzEV, szEV;
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
    4.1.3 18, 4.2.3 20,4.3.3 22/,
    zEVmax(s),
    xEVmax(j,s);
table a(i,j)
    1   2   3
1   30  28  27
2   22  23  21
3   31  34  34
4   24  20  24;
Variable z;
Positive Variables x(j);
Equations ucelfce, omez0(i), omez1(i), omez2(j), omez3(j);
ucelfce..           z =E= sum( j, c(j)*x(j) );
omez0(i)..          sum( j, a(i,j) * x(j) ) =L= b_UP(i);
omez1(i)..          sum( j, a(i,j) * x(j) ) =G= b_DOWN(i);
omez2(j)..          x(j) =L= x_UP(j);
omez3(j)..          x(j) =G= x_DOWN(j); #X.LOW(j) = X(j) zdrojak IS -> Zrovnat si tabulator, dvojstupnove prohlednut si tie smesi, ja za jaku dodoatecnu cenu umoznumujem posunutie Ipsilon
model vyroba / ucelfce, omez0, omez1, omez2, omez3 /;
solve vyroba maximizing z using LP;
display z.L, x.L;
*---------------------------------------------------------------------------
*  >>> neatly formatted output file  <<<
*  https://www.youtube.com/watch?v=erG5GmNWbv8
*---------------------------------------------------------------------------
file out / "vysledky.txt" /;
put out;
put "Vysledky a vstupy:" /;
put "==================" / /;

* ---------- column headings ----------
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

loop(j, put x_UP(j):16:2;);
loop(j, put x_DOWN(j):16:2;);
loop(i, loop(j,put a(i,j):15:2;););
put /;
put /;
#To Find the current value of X and cap the fucker
x.Lo(j) = x.L(j);
x.Up(j) = x.L(j);
loop(s,
    a(i,j) = as(i,j,s);
    
    solve vyroba maximizing z using LP;
    display z.L, x.L;
    
    zEVmax(s) = z.L;
    
    if ((vyroba.modelstat = 4) or (vyroba.modelstat = 19),
        zEVmax(s) = -INF
    );
    
    put s.tl:3, p(s):5:3, vyroba.modelstat:3:0, vyroba.solvestat:7:0, zEVmax(s):12:2;
    loop(j,
        put x.l(j):10:2;
        put "":2;
        put$(ord(j)=1) "":1;
        put$(ord(j)=2) "":1;
    );
    loop(i,
        put b_UP(i):13:2;
        put$(ord(i)=1) "":3;
        put$(ord(i)=2) "":3;
        put$(ord(i)=3) "":3;
        );
    loop(i, put b_DOWN(i):16:2;);
    loop(j, put x_UP(j):16:2;);
    loop(j, put x_DOWN(j):16:2;);
    loop(i, loop(j, put a(i,j):15:2;););
    put /;
);
put /;
* ---------- risk measures ----------
EzEV   = sum(s,p(s) * zEVmax(s) );
put "EzEV   = ", EzEV:12:2 /;
If(EzEV GT -INF,varzEV = sum(s,p(s) * zEVmax(s) * zEVmax(s) ) - EzEV * EzEV;);
If(EzEV EQ -INF,varzEV = UNDF;);
put "varzEV = ", varzEV:12:2 /;
szEV   = sqrt(varzEV);
put "stdzEV = ", szEV:12:2 /;