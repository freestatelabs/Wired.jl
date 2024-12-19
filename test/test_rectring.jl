using Wired
circ, rect = loadrings("testdata/testring.csv")

using Plots 

# Nrows = 100
# Ncols = 100
# xrange = range(0, 2, Nrows)
# zrange = range(-1, 1, Ncols)
# y = 0
# B = zeros(Nrows*Ncols,3)
# nodes = zeros(size(B))

# it = 1
# for i in range(1, Nrows)
#     for j in range(1, Ncols)
#         nodes[it,:] = [xrange[i], y, zrange[j]]
#         global it +=1
#     end
# end



circ = biotsavart!(B, nodes, rect)
Bmag = zeros(Nrows*Ncols)
Bmag = normrows!(Bmag, B)

function meshgrid(x, y, bmag)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    Bmag = []
    return X, Y
end

plot([x for x in xrange], B)


