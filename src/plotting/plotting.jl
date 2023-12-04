using Makie

function plotPoints!(fig::Figure, points::Vector{<:Vector{<:Real}})
    # Extract x and y coordinates from points
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]

    # Create an axis
    ax = Axis(fig[1, 1], aspect = DataAspect(), limits = (nothing, (-0.5, 1.5)), xticks= -5:0.5:5, xminorticksvisible = true )

    # Create a scatter plot
    scatter!(ax, x_coords, y_coords)

    # Return the figure
    fig
end

function plotPoints!(fig::Figure, ax::Axis, points::Vector{<:Vector{<:Real}})
    # Extract x and y coordinates from points
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]

    # Create a scatter plot
    scatter!(ax, x_coords, y_coords)

    # Return the figure
    fig
end

# Function to plot a Ray
function plotRay!(fig::Figure, ray::Ray, length::Real=1.0)
    # Calculate end point of the ray
    end_point = ray.origin .+ length .* (ray.direction ./ norm(ray.direction))

    # Create an axis with 1:1 aspect ratio
    ax = Axis(fig[1, 1], aspect = DataAspect(), limits = (nothing, (-0.5, 1.5)), xticks= -5:0.5:5, xminorticksvisible = true )

    # Plot the ray
    lines!(ax, [ray.origin[1], end_point[1]], [ray.origin[2], end_point[2]])

    # Return the figure
    fig
end

# Function to plot a LRay
function plotLRay!(fig::Figure, lray::LRay, length::Real=1.0)
    # Calculate end points of the LRay
    left_end_point = lray.origin .- length .* (lray.left ./ norm(lray.left))
    right_end_point = lray.origin .+ length .* (lray.right ./ norm(lray.right))

    # Create an axis with 1:1 aspect ratio
    ax = Axis(fig[1, 1], aspect = DataAspect(), limits = (nothing, (-0.5, 1.5)), xticks= -5:0.5:5, xminorticksvisible = true )

    # Plot the LRay
    lines!(ax, [left_end_point[1], lray.origin[1], right_end_point[1]], [left_end_point[2], lray.origin[2], right_end_point[2]])

    # Return the figure
    fig
end

# Function to plot an array of LRays
function plotLRays!(fig::Figure, lrays::Vector{LRay}, length::Real=1.0)
    # Create an axis with 1:1 aspect ratio
    ax = Axis(fig[1, 1], aspect = DataAspect(), limits = (nothing, (-0.1, 2)), xticks= -5:0.5:5, xminorticksvisible = true )

    # Loop over all LRays and plot them
    for lray in lrays
        # Calculate end points of the LRay
        left_end_point = lray.origin .- length .* (lray.left ./ norm(lray.left))
        right_end_point = lray.origin .+ length .* (lray.right ./ norm(lray.right))

        # Plot the LRay
        lines!(ax, [left_end_point[1], lray.origin[1], right_end_point[1]], [left_end_point[2], lray.origin[2], right_end_point[2]])
    end

    # Return the figure
    fig, ax
end
