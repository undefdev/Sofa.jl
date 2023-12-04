using Dates
using Serialization

Batches = Array{SubdividedSubBox,1}
BatchesBundle = Tuple{Batches, Int}

const TRIPLES :: TripleArray = [
    ( 697, 696, 985 ),
    (  7,   24,  25 ),
    (  33,  56,  65 ),
    (  56,  33,  65 ),
    (  24,   7,  25 ),
]

const GERVER = 2219531668871967//
               1000000000000000 # Gerver's constant, rounded down

function getCheckpointFiles(foldername)
    files = [f for f in readdir(foldername) if occursin(".dat", f)]
    # find most recent timestamps in filenames
    timestamps = [DateTime(Base.split(Base.split(f, ".")[1], "_")[end]) for f in files]
    # sort files by timestamp
    perm = sortperm(timestamps, rev=true)
    files = files[perm]
    return files
end

function load_or_initialize_state(foldername) :: BatchesBundle
    files = getCheckpointFiles(foldername)
    
    # If no batch files are found, return the initial state
    if isempty(files)
        return [[(subBoxFromTriples(TRIPLES), Rational{BigInt}(Inf), false)]], 1
    end

    batches, last_iteration = deserialize(joinpath(foldername, files[1]))
    return batches, last_iteration + 1
end

function processBatch( batch::SubdividedSubBox, triples::TripleArray, lowerBound::Real ) :: SubdividedSubBox
    # Process each box in the batch
    tasks = map( batch ) do s
        subdivideAndFilter( s, triples, lowerBound )
    end
    
    # Fetch results from each task and combine them
    refined = mapreduce( fetch, vcat, tasks )

    # Sort the refined subdivision in descending order based on the area
    sort!( refined, by = x -> x[2], rev = true )
    
    return refined
end

function processBatchWithBoxes( batch::Union{SubdividedSubBox,RefinementSubdividedSubBox}, triples::TripleArray, rBox::Vector{SubBox}, iTriples::TripleArray, lowerBound::Real ) :: RefinementSubdividedSubBox
    # Process each box in the batch
    tasks = map( batch ) do s
        subdivideWithRefinementAndFilter( s, triples, rBox::Vector{SubBox}, iTriples::TripleArray, lowerBound )
    end
    
    # Fetch results from each task and combine them
    refined = mapreduce( fetch, vcat, tasks )

    # Sort the refined subdivision in descending order based on the area
    sort!( refined, by = x -> x[2], rev = true )
    
    return refined
end

# function computeLowerBound( batch::SubdividedBox, triples::TripleArray, currentLowerBound::Rational ) :: Real
#     newLowerBound = currentLowerBound
#     for (box, _) in batch
#         newLowerBound = max( newLowerBound, centerAreaFromBoxAndTriples( box, triples ) )
#     end
#     return newLowerBound
# end

function computeBounds(foldername, lowerBound=GERVER, batchSize=2::Integer)
    if !isdir(foldername)
        mkdir(foldername)
    end

    batches, i = load_or_initialize_state(foldername)

    println( "writing to folder: ", foldername )
    # println( "computing lower bound..." )
    # for batch in batches
        # lowerBound = computeLowerBound( batch, TRIPLES, lowerBound )
    # end
    println( "lower bound: ", lowerBound )
    println( "batch size: ", batchSize )
    println( "### triples:")
    display( TRIPLES )

    activeIndex = findmax(x -> x[1][2], batches)[2]
    
    while true
        local sBox = processBatch(batches[activeIndex], TRIPLES, lowerBound)
        println("iteration :", i)
        println("   active index: ", activeIndex)

        # println("   recomputing lower bound..." )
        # lowerBound = computeLowerBound( sBox, TRIPLES, lowerBound )
        # println("   lower bound: ", lowerBound )
        # println("   lower bound approx: ", Float64(round(lowerBound, digits=8)))

        local currentBatchSize = min( batchSize, length( sBox ) )
        firstBatch, rest = sBox[1:currentBatchSize], sBox[currentBatchSize+1:end]
        if isempty(firstBatch)
            nBatches = length(batches)
            println("   batch is empty, reordering ", nBatches, " -> ", activeIndex )
            batches[activeIndex] = pop!(batches)
        else
            batches[activeIndex] = firstBatch

            if !isempty(rest) # if firstBatch is empty, then rest is empty too
                push!(batches, rest)
            end
        end

        # find index of largest area
        maxArea, activeIndex = findmax(x -> x[1][2], batches) # assumes batches are sorted in descending order

        println("   boxes in subdivision after iteration: ", sum(length.(batches)))
        println("   number of batches: ", length(batches))
        println("   currently largest area: ", maxArea)
        println("   currently largest area approx: ", Float64(round(maxArea, digits=8)))
        nextTarget = Float64(floor(maxArea, digits=2))
        # find indices of batches where largest elements have area greater or equal to the next lower second digit
        remainingBatches = findall(x -> x[1][2] >= nextTarget, batches)
        println("   batches greater than ", nextTarget, ": ", length(remainingBatches))

        # Every 1k iterations save batches to timestamped file
        if i % 1000 == 0
            println("#")
            println("# saving batches to file...")
            filename = "$(i)_$(now()).dat"
            filepath = joinpath(foldername, filename)
            serialize(filepath, (batches, i))
            # double check that serialization worked
            doublecheck, doublecheck_i = deserialize(filepath)
            @assert doublecheck_i == i
            maxArea_doublecheck, _ = findmax(x -> x[1][2], doublecheck)
            @assert maxArea == maxArea_doublecheck
            println("#")
            println("# saved batches to file: ", filename)
            println("#")
            # delete old files, keeping only the last 5
            files = getCheckpointFiles(foldername)
            # delete all but the 3 newest files
            for f in files[4:end]
                rm(joinpath(foldername, f))
            end
        end
        GC.safepoint()
        i += 1
    end
end

function refineBounds( box::SubBox, triples::TripleArray, lowerBound::Real, iterations=1::Integer, batchSize=2::Integer ) :: SubdividedSubBox
    sBox :: SubdividedSubBox = [(box,Rational{BigInt}(Inf),false)]
    batches = [sBox]
    local i=1
    activeIndex = 1
    while i!=iterations + 1 # allows us to use negative integers to indicate infinite iterations
        display("active index: $activeIndex")
        local sBox = processBatch( batches[activeIndex], triples, lowerBound)
        local currentBatchSize = min( batchSize, length( sBox ) )
        firstBatch, rest = sBox[1:currentBatchSize], sBox[currentBatchSize+1:end]
        batches[activeIndex] = firstBatch
        if !isempty( rest )
            push!( batches, rest )
        end
        # find index of largest area
        maxArea, activeIndex = findmax( x->x[1][2], batches )

        println("iteration:", i)
        println("   boxes in subdivision after iteration: ", sum(length.(batches)))
        println("   number of batches: ", length(batches))
        println("   currently largest area: ", maxArea)
        println("   currently largest area approx: ", Float64(round(maxArea, digits=8)))
        i+=1
        GC.safepoint()
    end
    return reduce(vcat, batches)
end

function refineBoundsWithBoxes( flat::Union{SubdividedSubBox,RefinementSubdividedSubBox}, triples::TripleArray, rBox::Vector{SubBox}, iTriples::TripleArray, lowerBound::Real, iterations=1::Integer, batchSize=2::Integer )# :: SubdividedSubBox
    # split flat into batches of size batchSize
    n = length( flat )
    batches :: Vector{Union{SubdividedSubBox,RefinementSubdividedSubBox}} = [ flat[i:min( i+batchSize-1, n )] for i in 1:batchSize:n ]
    local i=1
    activeIndex = 1
    while i!=iterations + 1 # allows us to use negative integers to indicate infinite iterations
        display("active index: $activeIndex")
        local sBox = processBatchWithBoxes( batches[activeIndex], triples, rBox::Vector{SubBox}, iTriples::TripleArray, lowerBound)
        local currentBatchSize = min( batchSize, length( sBox ) )
        firstBatch, rest = sBox[1:currentBatchSize], sBox[currentBatchSize+1:end]
        println("    improvement on batch ≈ ", Float64(round(batches[activeIndex][1][2], digits=8)) - Float64(round(firstBatch[1][2], digits=8)))
        batches[activeIndex] = firstBatch
        if !isempty( rest )
            push!( batches, rest )
        end
        # find index of largest area
        maxArea, activeIndex = findmax( x->x[1][2], batches )

        println("iteration:", i)
        println("   boxes in subdivision after iteration: ", sum(length.(batches)))
        println("   number of batches: ", length(batches))
        println("   currently largest area: ", maxArea)
        println("   currently largest area approx: ", Float64(round(maxArea, digits=8)))
        i+=1
        GC.safepoint()
    end
    return reduce(vcat, batches)
end

