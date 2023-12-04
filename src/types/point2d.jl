
Point2d = SVector{2, <:Real}
const HORIZONTAL = SVector{2,Rational{Int128}}(1, 0)

# rotates a point about the origin
rotate( p::Point2d, α::Real ) :: Point2d =
    [ p[1] * cos(α) - p[2] * sin(α), p[1] * sin(α) + p[2] * cos(α) ]

# rotates a point about the origin using a Pythagorean triple (a, b, c)
function rotate(p::Point2d, triple::Tuple{Int, Int, Int}) :: Point2d
    a, b, c = triple
    cosα = b // c
    sinα = a // c
    SA[ p[1] * cosα - p[2] * sinα, p[1] * sinα + p[2] * cosα ]
end

rotate90CW( p::Point2d ) :: Point2d = [ p[2], -p[1] ]
