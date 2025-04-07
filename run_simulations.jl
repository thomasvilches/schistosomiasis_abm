
#= 
run1(0)
run2(0)
run3(0)
 =#


# This is HTX because the wld is 2
#= 
run1(2,1,true,[0.0],[1],:dg,0.8,0)
run2(2,1,true,[0.0],[1],:dg,0.8,0)
run3(2,1,true,[0.0],[1],:dg,0.8,0)

run1(3,1,true,[0.0],[1],:dg,1.0,1)
run2(3,1,true,[0.0],[1],:dg,1.0,1)
run3(3,1,true,[0.0],[1],:dg,1.0,1)


run1(2,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run2(2,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run3(2,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)

run1(3,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run2(3,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run3(3,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
=#

#=

# This is KK because the wld is 2

run1(4,2,true,[0.0],[1],:dg,0.8,0)
run2(4,2,true,[0.0],[1],:dg,0.8,0)
run3(4,2,true,[0.0],[1],:dg,0.8,0)
run1(5,2,true,[0.0],[1],:dg,1.0,1)
run2(5,2,true,[0.0],[1],:dg,1.0,1)
run3(5,2,true,[0.0],[1],:dg,1.0,1)


run1(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run2(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)
run3(4,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0)

run1(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run2(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
run3(5,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1)
=#



# This is MDA
#= 
run1(6,2,true,[0.0],[1],:total,0.8,0)
run2(6,2,true,[0.0],[1],:total,0.8,0)
run3(6,2,true,[0.0],[1],:total,0.8,0)

run1(7,2,true,[0.0],[1],:total,1.0,1)
run2(7,2,true,[0.0],[1],:total,1.0,1)
run3(7,2,true,[0.0],[1],:total,1.0,1)


run1(6,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0)
run2(6,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0)
run3(6,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0)

run1(7,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1)
run2(7,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1)
run3(7,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1)
 =#


 for kk in 0.1:0.1:0.9

    println(string("kk ", kk))
    #= 
        run1(8; ks = true, pks = kk)
        run2(8; ks = true, pks = kk)
        run3(8; ks = true, pks = kk)
    =#
    #= 
    run1(9,1,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)
    run2(9,1,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)
    run3(9,1,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)

    run1(10,1,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)
    run2(10,1,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)
    run3(10,1,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)


    run1(9,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)
    run2(9,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)
    run3(9,1,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)

    run1(10,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
    run2(10,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
    run3(10,1,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
     =#

    run1(11,2,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)
    run2(11,2,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)
    run3(11,2,true,[0.0],[1],:dg,0.8,0; ks = true, pks = kk)

    run1(12,2,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)
    run2(12,2,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)
    run3(12,2,true,[0.0],[1],:dg,1.0,1; ks = true, pks = kk)


    run1(11,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)
    run2(11,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)
    run3(11,2,true,[0.5;1;2],[2;4;6;8;10],:dg,0.8,0; ks = true, pks = kk)

    run1(12,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
    run2(12,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
    run3(12,2,true,[0.5;1;2],[2;4;6;8;10],:dg,1.0,1; ks = true, pks = kk)
    

    run1(13,2,true,[0.0],[1],:total,0.8,0; ks = true, pks = kk)
    run2(13,2,true,[0.0],[1],:total,0.8,0; ks = true, pks = kk)
    run3(13,2,true,[0.0],[1],:total,0.8,0; ks = true, pks = kk)

    run1(14,2,true,[0.0],[1],:total,1.0,1; ks = true, pks = kk)
    run2(14,2,true,[0.0],[1],:total,1.0,1; ks = true, pks = kk)
    run3(14,2,true,[0.0],[1],:total,1.0,1; ks = true, pks = kk)


    run1(13,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0; ks = true, pks = kk)
    run2(13,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0; ks = true, pks = kk)
    run3(13,2,true,[0.5;1;2],[2;4;6;8;10],:total,0.8,0; ks = true, pks = kk)

    run1(14,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1; ks = true, pks = kk)
    run2(14,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1; ks = true, pks = kk)
    run3(14,2,true,[0.5;1;2],[2;4;6;8;10],:total,1.0,1; ks = true, pks = kk)
    
end

 
