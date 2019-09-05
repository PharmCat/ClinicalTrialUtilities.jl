struct RandomTable
    seq
    table
    t
end

function randomtable(;blocksize = 4, subject = 10, group = 2, ratio = [1,1], grseq = [["A", "B"], ["B", "A"]], seed = 1234)

    rng = MersenneTwister()
    if seed == 0  Random.seed!(rng) else Random.seed!(seed) end

    r = sum(ratio)

    if length(ratio) != group throw(ArgumentError("Length of ratio should be equlal group")) end
    if blocksize%r != 0 throw(ArgumentError("Sum of ratio should be in blocksize")) end
    if subject%r != 0 throw(ArgumentError("Sum of ratio should be in subject")) end

    if length(unique(length.(grseq))) != 1 throw(ArgumentError("Unequal grseq")) end

    rm = blocksize/r                                                            #ratio multiplication
    block = []                                                                  #block dummy
    for i = 1:group
        append!(block, fill(i, Int(rm*ratio[i])))
    end
    fbn = subject÷blocksize                                                     #block number
    rand  = []                                                                  #rand list
    for i = 1:fbn
        append!(rand, block[sample(1:blocksize, blocksize, replace=false)])     #generate and append random block
    end

    last = subject%blocksize                                                    #not full block
    if last != 0
        block = []
        rm = last/r                                                             #ratio multiplication
        for i = 1:group
            append!(block, fill(i, Int(rm*ratio[i])))
        end
        append!(rand, block[sample(1:last, last, replace=false)])
    end

    #seqrand = hcat(grseq[rand]...)

    return rand

end
