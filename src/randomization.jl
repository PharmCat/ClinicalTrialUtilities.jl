struct RandomTable
    seq
    table
    t
end
"""
    randomseq(;blocksize = 4, subject = 10, group = 2, ratio = [1,1],  seed = 0)

Randomization list generaton.

- ``blocksize`` - size of block;
- ``subject`` - subject number;
- ``group`` - number of groups;
- ``ratio`` - group ration, [1,1] means 1:1, [1,2] means 1:2, ets. ``length(ratio)`` should be equal ``group``;
- ``seed`` - RNG seed.
"""
function randomseq(;blocksize::Int = 4, subject = 10, group = 2, ratio = [1,1], seed = 0)

    #rng = MersenneTwister()
    if seed != 0 Random.seed!(seed) end

    r = sum(ratio)

    if length(ratio) != group throw(ArgumentError("Length of ratio should be equlal group")) end
    if blocksize%r != 0 throw(ArgumentError("Sum of ratio should be in blocksize")) end
    if subject%r != 0 throw(ArgumentError("Sum of ratio should be in subject")) end

    #if length(unique(length.(grseq))) != 1 throw(ArgumentError("Unequal grseq")) end

    rm = blocksize/r                                                            #ratio multiplication
    block = []                                                                  #block dummy
    for i = 1:group
        append!(block, fill(i, Int(rm*ratio[i])))
    end
    fbn = subject÷blocksize                                                     #block number
    rand  = Vector{Int}(undef, 0)                                                                  #rand list
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

"""
    randomtable(;blocksize = 4, subject = 10, group = 2,
        ratio = [1,1], grseq = ["AB", "BA"], seed = 0)


Randomization table generaton.

- ``blocksize`` - size of block;
- ``subject`` - subject number;
- ``group`` - number of groups;
- ``ratio`` - group ration, [1,1] means 1:1, [1,2] means 1:2, ets. ``length(ratio)`` should be equal ``group``;
- ``grseq`` - treatment sequence in each group;
- ``seed`` - RNG seed.
"""
function randomtable(;blocksize::Int = 4, subject::Int = 10, group::Int = 2, ratio::Vector = [1,1], grseq::Vector{String} = ["AB", "BA"], seed = 0)

    if length(unique(length.(grseq))) != 1 throw(ArgumentError("Unequal grseq")) end

    r = randomseq(;blocksize = blocksize, subject = subject, group = group, ratio = ratio,  seed = seed)
    subj = collect(1:length(r))
    seqrand  = vcat(grseq[r]...)
    nl = length(grseq[1]) + 3
    nm = Vector{Symbol}(undef, nl)
    np = Vector{Vector}(undef, nl)
    nm[1] = :Subject
    nm[2] = :Group
    nm[3] = :Sequence
    np[1] = subj
    np[2] = r
    np[3] = seqrand
    for i=1:length(grseq[1])
        nm[i+3] = Symbol("Period_"*string(i))
        np[i+3] = getindex.(seqrand, i)
    end
    return NamedTuple{Tuple(nm)}(Tuple(np))

    seqrand  = hcat(collect(1:length(r)), hcat(r,  seqrand))
    return seqrand
    df       = DataFrame(seqrand)
    rename!(df, [:Subject, :Group, :Sequence])
    for i=1:length(grseq[1])
        df[!, Symbol("Period_"*string(i))] = getindex.(df[!, :Sequence], i)
    end
    return df
end
