
struct Interval
    lo::Rational
    hi::Rational
end

"""
    A SubInterval is an Interval with an associated address, which is a BitVector.
    The length of the BitVector encodes how many times the original interval was split,
    while the value of each entry in the BitVector encodes whether the left or right subinterval was chosen.
    
    For example for a root Interval(0,1), the SubInterval(Interval(0,1), BitVector([])) is the original interval,
    while SubInterval(Interval(0,1), BitVector([0,1,0])) represents Interval(1//4,1//4 + 1//8).
"""
struct SubInterval
    root::Interval
    address::BitVector
    SubInterval( root::Interval, address=BitVector() ) = new( root, address )
end

function ⊆( i1 :: SubInterval, i2 :: SubInterval ) :: Bool
    n, m = length(i1.address), length(i2.address)
    return i1.root==i2.root && m<=n && i1.address[1:m] == i2.address
end

function boundsFromSubInterval(subinterval::SubInterval) :: Tuple{Rational, Rational}
    address = subinterval.address
    
    # Handle the case when the address is empty
    if isempty(address)
        return subinterval.root.lo, subinterval.root.hi
    end

    # Function to estimate bit length safely
    bit_length = x -> x == 0 ? 0 : ceil(Int, log2(abs(x)) + 1)
    
    # Estimate the maximum possible bit length required
    root = subinterval.root
    max_bit_length = max(bit_length(numerator(root.lo)), bit_length(denominator(root.lo)), 
                         bit_length(numerator(root.hi)), bit_length(denominator(root.hi))) + length(address)
    
    # Decide whether to use BigInt based on the estimated bit length
    use_bigint = max_bit_length > 62
    T = use_bigint ? (max_bit_length > 126 ? BigInt : Int128) : Int64
    den = one(T) << length(address)

    offset = sum( T(bit) << (length(address) - i) for (i, bit) in enumerate(address) ) // den

    # Now scale and translate the result to the root interval
    root_size = root.hi - root.lo
    lower_bound = root.lo + root_size * offset
    upper_bound = lower_bound + root_size // den

    return lower_bound, upper_bound
end

mid( subInterval::SubInterval ) :: Rational = sum( boundsFromSubInterval( subInterval ) )//2

lerp( interval::Interval, t::Rational ) :: Rational = interval.lo + (interval.hi - interval.lo) * t

function split( subinterval::SubInterval ) :: NTuple{2,SubInterval}
    root, address = subinterval.root, subinterval.address
    return SubInterval( root, vcat( address, false ) ),
           SubInterval( root, vcat( address, true ) )
end
