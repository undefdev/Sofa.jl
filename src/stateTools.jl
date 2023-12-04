

function batchSubdivision( subdivision::SubdividedSubBox, batchSize=2::Integer ) :: Batches
    n = length( subdivision )
    return [subdivision[i:min(i+batchSize-1, n)] for i in 1:batchSize:n]
end

function splitStateFileAtCutoff( filePath::String, cutoff::Float64, batchSize=2::Integer ) :: Tuple{Batches, Batches, Int}
    # Load the state file
    batches, iteration = deserialize( filePath )

    # Flatten the batches
    flatData = reduce( vcat, batches )

    # Sort the flattened data based on the area
    sort!( flatData, by = x -> x[2], rev = true )

    # Find the index where the area falls below the cutoff
    cutoffIndex = findfirst( x -> x[2] < cutoff, flatData )

    # Handle case where all elements are above the cutoff
    if cutoffIndex === nothing
        return flatData, [], iteration
    end

    # Split the data at the cutoff index
    greaterEqCutoff = flatData[1:cutoffIndex - 1]
    lessThanCutoff = flatData[cutoffIndex:end]

    # Re-batch the split data
    greaterEqBatches = batchSubdivision( greaterEqCutoff, batchSize )
    lessBatches      = batchSubdivision( lessThanCutoff, batchSize )

    return greaterEqBatches, lessBatches, iteration
end

function splitAndBatch( flatData::SubdividedSubBox, n_parts::Int, batch_size::Integer ) :: Array{Batches, 1}
    # Initialize an array of empty arrays for each part
    parts = [Vector{typeof(flatData[1])}() for _ in 1:n_parts]

    # Distribute the flattened data across the parts
    for (i, item) in enumerate( flatData )
        push!( parts[(i - 1) % n_parts + 1], item )
    end

    # Batch each part and return
    return batchSubdivision.( parts, batch_size )
end

function isEqualSubdivision( subdivision1::SubdividedSubBox, subdivision2::SubdividedSubBox ) :: Bool
    # Check if the lengths are equal
    if !(length( subdivision1 ) == length( subdivision2 ))
        println("Length mismatch")
        return false
    end
    
    subBoxes, loadedSubBoxes = (x[1] for x in subdivision1), (x[1] for x in subdivision2)
    areas, loadedAreas       = (x[2] for x in subdivision1), (x[2] for x in subdivision2)

    if any( a1!=a2 for (a1, a2) in zip( areas, loadedAreas ) )
        println("Area mismatch")
        return false
    end
    
    for (s1, s2) in zip( subBoxes, loadedSubBoxes )
        if any( i1.root!=i2.root || i1.address!=i2.address for (i1, i2) in zip( s1, s2 ) )
            println("SubInterval mismatch")
            return false
        end
    end

    # If we reach here, then the subdivisions are equal
    return true
end
function serializeAndVerifyBatches( batches::Batches, iteration::Int, directory::String )
    # Ensure the directory exists
    isdir( directory ) || mkdir( directory )

    # Construct the file path
    filePath = joinpath( directory, "$(iteration)_$(Dates.now()).dat" )

    # Serialize the batch along with iteration number
    serialize( filePath, (batches, iteration) )

    # Deserialize to verify
    loadedBatches, loadedIteration = deserialize( filePath )

    if !(loadedIteration == iteration)
        error( "Iteration number mismatch for batches at $filePath" )
    end

    verified = all( (isEqualSubdivision( b1, b2 ) for (b1, b2) in zip( batches, loadedBatches )) )
    if !verified
        error( "Subdivision mismatch for batches at $filePath" )
    end

    println( "All batches successfully serialized and verified." )
end

function mergeBatches( batchesArray::Array{Batches, 1} ) :: SubdividedSubBox
    return sort!( reduce( vcat, reduce( vcat, batchesArray ) ), by=x->x[2], rev=true )
end

"""
takes a foldername and a batch size and splits all files in the folder into batches of the given size
folder should contain subfolders `big` and `small` which contain the files to be split
files in `small` don't need to be split, so they are just copied over
"""
function resplitFiles( foldername::String, cutoff::Float64, nParts::Integer, batchSize=2::Integer, lastIteration=0::Integer )
    # Get all files in the folders
    bigFiles   = readdir( joinpath( foldername, "big" ) )
    smallFiles = readdir( joinpath( foldername, "small" ) )

    # Filter out files that don't end with .dat
    bigFiles   = filter( x -> endswith( x, ".dat" ), bigFiles )
    smallFiles = filter( x -> endswith( x, ".dat" ), smallFiles )

    # Filter out files that don't start with a number
    bigFiles   = filter( x -> isdigit( x[1] ), bigFiles )
    smallFiles = filter( x -> isdigit( x[1] ), smallFiles )

    # Split the files at the cutoff
    greaterEqCutoffs, lessThanCutoffs = [], []
    for (greater, less, i) in (splitStateFileAtCutoff( joinpath( foldername, "big", bigFile ), cutoff, batchSize ) for bigFile in bigFiles)
        push!( greaterEqCutoffs, (greater, i - lastIteration) )
        push!( lessThanCutoffs,  (less,    i - lastIteration) )
    end

    for (small, i) in (deserialize( joinpath( foldername, "small", smallFile ) ) for smallFile in smallFiles)
        push!( lessThanCutoffs, (small, i - lastIteration) )
    end

    bigIteration   = lastIteration + sum( x[2] for x in greaterEqCutoffs )
    smallIteration = lastIteration + sum( x[2] for x in lessThanCutoffs )

    bigBatches = mergeBatches( [x[1] for x in greaterEqCutoffs] )
    smallBatches = mergeBatches( [x[1] for x in lessThanCutoffs] )

    # Re-batch the split files
    bigBatchesArray = splitAndBatch( bigBatches, nParts, batchSize )
    smallBatchesArray = [smallBatches]

    outputFolderName = joinpath( foldername, "resplit_$(Dates.now())" )
    if !isdir( outputFolderName )
        mkdir( outputFolderName )
        mkdir( joinpath( outputFolderName, "small" ) )
    end
    # Serialize the batches
    for (i, bs) in enumerate( bigBatchesArray )
        dir = joinpath( outputFolderName, "big_$(i)" )
        if !isdir( dir )
            mkdir( dir )
        end
        serializeAndVerifyBatches( bs, bigIteration, dir )
    end
    dir = joinpath( outputFolderName, "small" )
    serializeAndVerifyBatches( smallBatchesArray, smallIteration, dir )
end