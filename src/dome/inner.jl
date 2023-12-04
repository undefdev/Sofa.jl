
function cornerIsContainedInOtherTriangle( ts::Array{Rational,4}, i::Integer, j::Integer ) :: Bool
    local ijlr, ijrl = ts[i,j,LEFT,RIGHT], ts[i,j,RIGHT,LEFT]
    return ijrl < 0 < ijlr # i is strictly contained within legs of j
end

function findNonContainedCorners( ts :: Array{Rational,4} ) :: Vector{Bool}
    n = size( ts, 1 )
    # find all corners that are not contained in another triangle
    nonContained =[ !cornerIsContainedInOtherTriangle( ts, i, j ) for i in 1:n, j in 1:n ]
    return [ reduce( &, nonContained, dims=2 )... ]
end

"""
    innerDome( innerLRays::SubArray{LRay, 1}, innerTs::Array{Rational, 4} ) :: Tuple{Array{Matrix{Rational},1}, Array{Matrix{Integer},1}}

    Compute the inner dome of the corridor intersection.
    The inner dome is the union of the triangles formed by the lower left and right legs of the corridors along with the x axis.
    
    The function returns the connected components of the inner dome, where each component is a matrix of inner dome vertices, as well as the indices of the LRays that form the inner dome.
"""
function innerDome( innerLRays::SubArray{LRay, 1}, innerTs::Array{Rational, 4} ) :: Tuple{Array{Matrix{Rational},1}, Array{Matrix{Integer},1}}
    # find all corners that are not contained in another triangle (peaks of the inner dome)
    peakMask = findNonContainedCorners( innerTs ) .& [ lRay.origin[2] > 0 for lRay in innerLRays ] # peaks are above x axis
    peakIndices = findall( peakMask )

    # if there are no peaks, there is no inner dome
    if isempty( peakIndices ) return [], [] end

    # sort peak indices by x coordinate
    sort!( peakIndices, by=peakIndex -> innerLRays[peakIndex].origin[1] )

    # compute intersection times for x-axis
    xAxisIntersectionTs = getIntersectionTsWithXAxis( innerLRays )

    # descend all the peaks left and right until they meet or reach the x axis, 0 indices delimit connected components
    descendsMatrix = computeDescendsMatrix( peakIndices, innerLRays, innerTs, xAxisIntersectionTs )

    # We now create a list of "index-components" [ component1, component2, ...],
    # where each component is a matrix where the first column is the indices of the LRays that form the component,
    # and the second column indicates whether the index is from the left or right leg of the corresponding LRay using LEFT and RIGHT.
    indexComponents = segment( descendsMatrix, xAxisIntersectionTs )
    
    # compute the vertices of each component and return them as list of matrices, along with the list of the corresponding "index-components"
    #components = verticesFromIndexComponent.( indexComponents, innerLRays, innerTs )
    components = [ verticesFromIndexComponent( indexComponent, innerLRays, innerTs ) for indexComponent in indexComponents ]

    return components, indexComponents
end

function innerDomeTruncated( innerLRays::SubArray{LRay, 1}, innerTs::Array{Rational, 4} ) :: Tuple{Array{Matrix{Rational},1}, Array{Matrix{Integer},1}} 
    components, indexComponents = innerDome( innerLRays, innerTs )
    return truncate.( components ), indexComponents
end

"""
    computeDescendsMatrix( peakIndices, innerLRays, innerTs, xAxisIntersectionTs ) :: Matrix{Integer}

    Descend all the peaks left and right until they meet or reach the x axis, 0 indices delimit connected components.
    Returns an n×2 matrix, where n is not known in advance.

    The left column contains the indices of the LRays that form the component, with 0 representing the x axis.
    The right column indicates whether the index is from the left or right leg of the corresponding LRay using LEFT and RIGHT.

    Example output:
    n×2 Matrix{Int64}:
     0  1
     2  1
     3  1
     ⋮  
     2  2
     3  2
     0  2
"""
function computeDescendsMatrix( peakIndices, innerLRays, innerTs, xAxisIntersectionTs::Array{Rational,3} ) :: Matrix{Integer}
    # Create an iterator that yields the submatrices
    descends_iter = ( let descendPartL = descend(innerLRays, innerTs, peakIndex, peakIndex, (@view xAxisIntersectionTs[:,2,LEFT]),  LEFT,  SLOW),
                          descendPartR = descend(innerLRays, innerTs, peakIndex, peakIndex, (@view xAxisIntersectionTs[:,2,RIGHT]), RIGHT, SLOW)
                        [ reverse!( descendPartL )  fill( LEFT,  length( descendPartL ) );
                                    descendPartR    fill( RIGHT, length( descendPartR ) ) ]
                      end for peakIndex in peakIndices )
      
    # Use reduce to concatenate all the submatrices into the descends matrix
    return reduce(vcat, descends_iter)
end

"""
    segment( descendsMatrix::Matrix{Integer} ) :: Array{SubArray{Integer,2},1}

    Assumes descendsMatrix contains at least one nonempty submatrix delimited by rows with a 0 in the first column.
    Returns these nonempty submatrices as an array of subarrays.
"""
function segment( descendsMatrix::Matrix{Integer}, xAxisIntersectionTs::Array{Rational,3} ) :: Array{SubArray{Integer,2},1}
    # find indices of rows with 0 in first column
    zeroIndices = findall( descendsMatrix[:,1] .== 0 )

    # get the nonempty intervals of indices between the 0s
    delimitedIntervals = [ let iStart = iZero+1,
                               iEnd   = iNextZero-1
                           iStart:iEnd 
                           end for (iZero, iNextZero) in neighbors( zeroIndices ) if iNextZero-iZero>1 ]
    

    # loop over intervals and merge them if they overlap
    # for this we first collect the indices of the intervals that we want to merge
    componentIntervals = [ [ delimitedIntervals[1] ] ]
    for interval in delimitedIntervals[2:end]
        lastInterval = componentIntervals[end][end]
        lastEndIndex = last( lastInterval )
        newStartIndex = first( interval )
        lastMaxXAxisIntersection    = xAxisIntersectionTs[descendsMatrix[lastEndIndex, 1],1,RIGHT]
        currentMinXAxisIntersection = xAxisIntersectionTs[descendsMatrix[newStartIndex,1],1,LEFT]

        # if intervals overlap, merge them by replacing the view in components[end] with the merged view            
        if lastMaxXAxisIntersection > currentMinXAxisIntersection
            push!( componentIntervals[end], interval  )
        else
            push!( componentIntervals, [ interval ] )
        end
    end

    # reduce each array of intervals to a list of indices, we use Iterators.map to avoid allocating a new array
    mergedComponentIndices = Iterators.map( intervals -> reduce( vcat, intervals ), componentIntervals )

    # map the indices to the corresponding submatrices
    return [ @view descendsMatrix[indices,:] for indices in mergedComponentIndices ]
end

function verticesFromIndexComponent( indexComponent::SubArray{Integer,2}, lRays::SubArray{LRay,1}, ts::Array{Rational,4} ) :: Matrix{<:Rational}
    n = size( indexComponent, 1 ) + 1 # +1 for the intersection with the x axis
    vertices = zeros( Rational, 2, n )

    # get the intersection points with the x axis
    vertices[1,1]   = getIntersectionTsWithXAxis( leftRay(  lRays[indexComponent[1,1]] ) )[1]
    vertices[1,end] = getIntersectionTsWithXAxis( rightRay( lRays[indexComponent[end,1]] ) )[1]

    # get the remaining intersection points
    for (idx, ((i,iSide),(j,jSide))) in enumerate( neighbors( eachrow( indexComponent ) ) )
        vertices[:,idx+1] = stepLRayByT( lRays[i], iSide, ts[i,j,iSide,jSide] )
    end

    return vertices
end
 