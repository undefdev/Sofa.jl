"""
    selectBox( mat::Array{SubInterval,2}, n::Int ) :: SubBox 

    Returns a `SubBox` that is the result of selecting the first or second `SubInterval` of each column of the given matrix, depending on the corresponding bit of the given integer.
"""
function selectBox( mat::Array{SubInterval,2}, n::Int )  :: SubBox
    k = size(mat, 2)
    # Initialize box as array of k+1 intervals, but don't initialize the intervals
    box = SubBox(undef, k+1)
    # set the first interval 0 and the rest to the corresponding interval in the matrix
    box[1] = SubInterval(Interval(0,0), BitVector())
    for i in 1:k
        box[i+1] = n & (1 << (i-1)) == 0 ? mat[1,i] : mat[2,i]
    end
    return box
end

function initializeSubdivision( box::SubBox, triples::TripleArray ) :: SubdividedSubBox
    k = length( box ) - 1 # number of intervals in the box (excluding the first interval, which is always [0,0])
    # Split every non-zero interval in the box and put the two intervals in a column of a matrix
    intervals = Matrix{SubInterval}(undef, 2, k)
    for i in 1:k
        intervals[:,i] .= split( box[i+1] ) # skipping first interval, since it's always [0,0]
    end

    subdivision = SubdividedSubBox(undef, 2^k)
    Threads.@threads for i in 0:2^k-1
        local box = selectBox( intervals, i )
        local area = minkowskiAreaFromBoxAndTriples( box, triples )
        subdivision[i+1] = (box, area, false)
    end
    return subdivision
end

function getLargestSubBoxes( subdivision )
    # If there is no subbox with area larger than the lower bound that hasn't already been handled, return an empty subdivision
    if isempty( subdivision ) return subdivision  end

    # Find maximal area of the subboxes
    maxArea = maximum( x->x[2], subdivision )

    # Return new subdivision that only contains boxes with area maxArea
    return filter( x->x[2]==maxArea, subdivision )
end
function filterSubdivision( subdivision::Vector{T}, lowerBound::Real ) :: Vector{T} where T
    sub = filter( x->lowerBound<=x[2], subdivision )
    return getLargestSubBoxes( sub )
end
function filterSubdivision( subdivision::Vector{T}, lowerBound::Real, upperBound::Rational ) :: Vector{T} where T
    sub = filter( x->lowerBound<=x[2]<upperBound, subdivision )
    return getLargestSubBoxes( sub )
end
function finalizeRefinementSubdivision( subdivision::RefinementSubdividedSubBox, box::SubBox )
    # # If there is no subbox with area larger than the lower bound that hasn't already been handled, return an empty subdivision
    if isempty( subdivision ) return subdivision  end

    # Only subboxes with maximal area remain, so we can choose the area of the first one
    maxArea = subdivision[1][2]

    # Push initial box with new upper bound to the subdivision, such that the subdivision can be further refined
    push!( subdivision, (box, maxArea, true, true) )

    return subdivision
end
"""
    subdivideAndFilter( box::SubBox ) :: SubdividedSubBox

    Returns a list of boxes that are the result of
        - splitting the given box at its midpoint into 2^k subboxes, where k is the number of intervals in the box
        - filtering out boxes that
            - have an area smaller than the lower bound
            - have already been extracted from this box
            - have an area smaller than the largest area of these subboxes

    Whether a subbox has already been extracted is determined by the box's upper bound, which is
        - `minkowskiAreaFromBoxAndTriples( box, triples )`      if the box hasn't been subdivided yet
        - the largest minkowski area of the subboxes            if the box has been subdivided
    """
function subdivideAndFilter( sub::SubBoxNode, triples::TripleArray, lowerBound::Real ) :: SubdividedSubBox
    box, upperBound, seen = sub # upperBound and seen are used to filter subboxes that have already been extracted from the subdivision

    # split box into 2^(k-1) subboxes (we don't need to split the first interval, since it's always [0,0])
    subdivision = initializeSubdivision( box, triples )

    # Filter out boxes that don't surpass the lower bound
    if seen
        # If this is not the first time we subdivide the box, we know that every box with area>=upperBound has already been extracted
        # Filter out boxes that have already been extracted from the subdivision
        subdivision = filterSubdivision( subdivision, lowerBound, upperBound )
    else
        subdivision = filterSubdivision( subdivision, lowerBound )
    end

    # # If there is no subbox with area larger than the lower bound that hasn't already been handled, return an empty subdivision
    if isempty( subdivision ) return subdivision  end

    # Only subboxes with maximal area remain, so we can choose the area of the first one
    maxArea = subdivision[1][2]

    # Push initial box with new upper bound to the subdivision, such that the subdivision can be further refined
    push!( subdivision, (box, maxArea, true) )

    return subdivision
end

function maxIntersectionAreaWithRefinementBox( box::SubBox, triples::TripleArray, rBox::Vector{SubBox}, rTriples::TripleArray ) :: Rational
    mergedTriples = [triples; rTriples] 
    return maximum( minkowskiAreaFromBoxAndTriples( [box; b], mergedTriples) for b in rBox )
end
function refinementSubdivisionFromSubdivision( subdivision::SubdividedSubBox, triples::TripleArray, rBox::Vector{SubBox}, rTriples::TripleArray ) :: RefinementSubdividedSubBox
    nSubboxes = length( subdivision )
    refinedSubdivision = RefinementSubdividedSubBox(undef, nSubboxes)
    
    Threads.@threads for i in 1:nSubboxes
        # compute refined area
        local box, _, _ = subdivision[i]
        local area = maxIntersectionAreaWithRefinementBox( box, triples, rBox, rTriples )
        refinedSubdivision[i] = (box, area, false, true)
    end

    return refinedSubdivision
end

function subdivideWithRefinementAndFilter( sub::SubBoxNode, triples::TripleArray, rBox::Vector{SubBox}, rTriples::TripleArray, lowerBound::Real ) :: RefinementSubdividedSubBox
    box, upperBound, seen = sub
    
    # if the box has not been subdivided, just convert it to a refinement box and return it
    if !seen
        # get the greatest area that can be obtained by intersecting with corridors represented by boxes in rBox
        local area = maxIntersectionAreaWithRefinementBox( box, triples, rBox, rTriples )
        return [(box, area, false, true)]
    end

    # oterwise we need to find out what boxes have already been extracted
    # and patch the seen box to have the correct bound
    
    # split box into 2^(k-1) subboxes (we don't need to split the first interval, since it's always [0,0])
    subdivision = initializeSubdivision( box, triples )
    
    # Filter out boxes that don't surpass the lower bound
    subdivision = filter( x->lowerBound<=x[2]<upperBound, subdivision )

    refinedSubdivision = refinementSubdivisionFromSubdivision( subdivision, triples, rBox, rTriples )
    
    refinedSubdivision = filterSubdivision( refinedSubdivision, lowerBound )
    
    return finalizeRefinementSubdivision( refinedSubdivision, box )
end
function subdivideWithRefinementAndFilter( sub::RefinementSubBoxNode, triples::TripleArray, rBox::Vector{SubBox}, rTriples::TripleArray, lowerBound::Real ) :: RefinementSubdividedSubBox
    box, upperBound, seen, rSeen = sub # upperBound and seen are used to filter subboxes that have already been extracted from the subdivision
    @assert rSeen

    # split box into 2^(k-1) subboxes (we don't need to split the first interval, since it's always [0,0])
    subdivision = initializeSubdivision( box, triples )

    # Filter out boxes that have already been extracted from the subdivision
    if seen
        subdivision = filter( x->lowerBound<=x[2]<upperBound, subdivision )
    else
        subdivision = filter( x->lowerBound<=x[2], subdivision )
    end

    refinedSubdivision = refinementSubdivisionFromSubdivision( subdivision, triples, rBox, rTriples )
    
    # Filter out boxes that don't surpass the lower bound
    if seen
        # If this is not the first time we subdivide the box, we know that every box with area>=upperBound has already been extracted
        # Filter out boxes that have already been extracted from the subdivision
        refinedSubdivision = filterSubdivision( refinedSubdivision, lowerBound, upperBound )
    else
        refinedSubdivision = filterSubdivision( refinedSubdivision, lowerBound )
    end
    
    return finalizeRefinementSubdivision( refinedSubdivision, box )
end

