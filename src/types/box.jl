# Box is an alias of an array of intervals
Box = Vector{Interval}
# SubBox is an alias of an array of subintervals
SubBox = Vector{SubInterval}

# SubBoxNode is a tuple of a SubBox, a rational number and a boolean value.
# The boolean value indicates whether the box has already been subdivided.
# The rational number is used to filter out subboxes that have already been extracted from the subdivision,
# if the box has already been subdivided.
# Initially, this rational number represents the largest area of a connected component
# of the MinkowskiCorridor associated with the box.
# As the box undergoes further subdivisions, this rational number is updated to serve
# as an adaptive upper threshold. This threshold is used to filter out sub-boxes that have
# been previously processed or explored, helping to narrow down the search space in the subdivision process.
# Importantly, the rational number always denotes an area of a connected component of the
# MinkowskiCorridor linked with a box that is contained within the paired box of that specific rational number.
SubBoxNode = Tuple{SubBox,Rational{BigInt},Bool} 
RefinementSubBoxNode = Tuple{SubBox,Rational{BigInt},Bool,Bool} 

# `SubdividedSubBox` is defined as an array of `SubBoxNode`s.
SubdividedSubBox = Vector{SubBoxNode}
RefinementSubdividedSubBox = Vector{RefinementSubBoxNode}

# function that splits a box into two boxes at the midpoint of the i-th interval
function split(box::SubBox, i::Int)
    box1 = copy(box)
    box2 = copy(box)
    box1[i], box2[i] = split( box[i] )
    return box1, box2
end

⊆( box1::SubBox, box2::SubBox ) :: Bool = all( x->x[1]⊆x[2], zip( box1, box2 ) )

function initialBoxFromTriples( triples::TripleArray ) :: Box 
    local n = length( triples )
    a₁, b₁, c₁ = triples[1]
    
    # initialize box as array of 2n intervals, but don't initialize the intervals
    box = Box(undef, 2n)
    
    sec₁ = c₁ // b₁
    box[2] = Interval( 0, sec₁ ) # lower bound could even be a bit higher than 0, as can be verified by computation 
    box[1] = Interval( 0, 0 )

    # get initial x₁, x₂ from the chosen triple index j=1
    csc₁ = c₁ // a₁
    x₁, x₂ = -csc₁ * (sec₁ + 1), sec₁

    for i in 2:n
        local I, J = 2i - 1, 2i
        a, b, c = triples[i]
        
        cosᵢ = b // c
        sinᵢ = a // c
        
        box[I] = Interval(
            x₁ * cosᵢ + (1 - cosᵢ) * cosᵢ // sinᵢ,
            x₂ * cosᵢ + sinᵢ - 1
        )
        box[J] = Interval(
            -x₂ * sinᵢ + (1 - sinᵢ) * sinᵢ // cosᵢ,
            -x₁ * sinᵢ + cosᵢ - 1
        )
    end

    return box
end

function subBoxFromTriples( triples::TripleArray ) :: SubBox
    box = initialBoxFromTriples( triples )
    return SubInterval.( box )
end

"""
    refinementBoxFromTriples( triples::TripleArray, nSplits::Integer, nDiscardTriples::Integer ) :: SubdividedSubBox

    Returns subboxes that are the result of splitting the initial box at the midpoint of each interval nSplits times.
    Skips the first nDiscardTriples triples.
"""
function refinementBoxesFromTriples( triples::TripleArray, nBisections::Integer, nDiscardTriples::Integer, lowerBound::Real ) :: Vector{SubBox}
    box = initialBoxFromTriples( triples )
    prefix = SubInterval.(box[1:2*nDiscardTriples])
    intervalsToSplit = box[2*nDiscardTriples+1:end]
    addresses = [BitVector(digits(i, base=2, pad=nBisections)) for i in 0:(2^nBisections - 1)]
    boxes = [ (x-> SubInterval( x[1], x[2] )).( comb ) for comb in combineElements( intervalsToSplit, addresses ) ]
    subboxes = []
    for b in boxes
        pbox = [prefix; b]
        area = minkowskiAreaFromBoxAndTriples( pbox, triples )
        if area > lowerBound
            push!( subboxes, b ) # discarding the prefix again
        end 
    end
    return subboxes
end
