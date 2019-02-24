
export expandPolygon, insidePolygon


function expandPolygon( x::Vector{Float64}, y::Vector{Float64},  # polygon points
	                     LL )  # distance
# Expand the polygon.  Each line segment goes out a distance of LL
# in the normal direction.
	
npts = length(x)

if length(y) != npts
   error("length(y) != npts")
end

xs = Array{Float64}(npts*2)
ys = Array{Float64}(npts*2)

idx = 1

for i = 1:npts-1
   
   theta = atan2(y[i+1] - y[i], x[i+1] - x[i])
   theta += pi/2  # rotate by 90 for normal direction.

   lct = LL * cos(theta)
   lst = LL * sin(theta) 

   xs[idx] = x[i] + lct
   ys[idx] = y[i] + lst 
   
   xs[idx+1] = x[i+1] + lct
   ys[idx+1] = y[i+1] + lst   

   idx += 2
end # i

# end point
theta = atan2(y[1] - y[npts], x[1] - x[npts])
theta += pi/2

lct = LL * cos(theta)
lst = LL * sin(theta) 

xs[idx] = x[npts] + lct
ys[idx] = y[npts] + lst   
   
xs[idx+1] = x[1] + lct
ys[idx+1] = y[1] + lst   
	
return xs, ys	
end # function expandPolygon

#---------------------------------------------------------------------

function insidePolygon( x::Vector{Float64}, y::Vector{Float64},  # polygon points
  	                     xp,yp)
#  Is point (xp,yp) inside polygon defined by x,y ?
	
n = length(x)

if length(y) != n
   error("length(y) != n")
end
	
inside = false	
	
for i = 1:n
   i1 = i + 1
   if i1 > n
      i1 = 1
   end
   dy = y[i1] - y[i]

   if dy != 0.0
      y1 = (yp <= y[i])
      y2 = (yp >  y[i1])

      if y1 == y2 

         rslope = (x[i1] - x[i]) / dy
         d = (xp-x[i]) - (yp-y[i]) * rslope
         
         if abs(d) < 1.e-4
         	 d = 0.0
         end

         if d <= 0.0
            inside = !inside
         end

      end  # y1 == y2 

   end # (dy != 0.0)

end # i
	
	
return inside	
end # function insidePolygon
