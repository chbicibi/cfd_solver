#!ruby -Ku

def naca4 mr, pr, t, c
  n = 10000
  m = mr * c
  p = pr * c
  ct5 = 5 * c * t
  yc = -> x do
    x < p ?
      m * x       / p       ** 2 * ( 2 * p -     x / c) :
      m * (c - x) / (1 - p) ** 2 * (-2 * p + 1 + x / c)
  end
  dyc_dx = -> x do
    x < p ?
      2 * m / p       ** 2 * (p - x / c) :
      2 * m / (1 - p) ** 2 * (p - x / c)
  end
  yt = -> x do
    xc = x / c
    ct5 * (  0.2969 * xc ** 0.5 +
            -0.1260 * xc        +
            -0.3516 * xc ** 2   +
             0.2843 * xc ** 3   +
            -0.1015 * xc ** 4 )
  end
  shape = -> x do
    ytx = yt[x]
    ycx = yc[x]
    tt = Math.atan dyc_dx[x]
    yts = ytx * Math.sin(tt)
    ytc = ytx * Math.cos(tt)
    [ x,
      ycx,
      x - yts,
      x + yts,
      ycx + ytc,
      ycx - ytc ]
  end
  n.times.map { |x| shape[x / (n - 1.0)] }
end

def rotate x, y, cx, cy, a
  [
    (x - cx) * Math.cos(a) - (y - cy) * Math.sin(a) + cx,
    (x - cx) * Math.sin(a) + (y - cy) * Math.cos(a) + cy
  ]
end

def grid_wing nx, ny, param, alpha
  a = ny.times.map { |iy| Array.new(nx, 1) }
  f = -> x, y { a[y.to_i][x.to_i] = 0 }
  scale = 200
  a_rad = (-Math::PI / 180) * alpha
  cx = nx / 10
  cy = ny / 2
  naca4(*param).map do |x, yc, xu, xl, yu, yl|
    ( scale*yu).to_i.times { |iy| f[*rotate(cx+scale*xu, cy-iy, cx + scale / 2, cy, a_rad)] }
    (-scale*yl).to_i.times { |iy| f[*rotate(cx+scale*xl, cy+iy, cx + scale / 2, cy, a_rad)] }
    # [x, yc, xu, xl, yu, yl, t]
  end
  a
end

def grid_circle nx, ny, cx, cy, r
  r2 = r ** 2
  ny.times.map do |iy|
    nx.times.map do |ix|
      (ix - cx) ** 2 + (iy - cy) ** 2 < r2 ? 0 : 1
    end
  end
end

def grid_arc nx, ny, r
  cx = nx / 2
  cy = Math.sqrt(r ** 2 - cx ** 2)
  grid_circle nx, ny, cx, -cy, r
end

def grid_plate nx, ny, t
  ny.times.map do |iy|
    nx.times.map do |ix|
      (ix > 0 and iy < t) ? 0 : 1
    end
  end
end

def grid_rect nx, ny, t, b, l, r
  ny.times.map do |iy|
    nx.times.map do |ix|
      (ix >= l and ix < r and iy < t and iy >= b) ? 0 : 1
    end
  end
end

def encode grid
  grid.map { |r| r.map(&:to_s).join(",") }
end

def main
  nx, ny = File.readlines("in2d.txt")[24].scan(/\d+/).map(&:to_i)
  # p grid_wing nx, ny
  # data = encode(grid_wing nx, ny, [0.02, 0.4, 0.12, 1.0])
  data = encode(grid_wing nx, ny, [0.0, 0.0, 0.12, 1.0], 0)
  # data = encode(grid_circle nx, ny, nx / 4, ny / 2, ny / 8)
  # data = encode(grid_arc nx, ny, nx * 2)
  # data = encode(grid_plate nx, ny, 1)
  # data = encode(grid_rect nx, ny, ny / 2, 0, ny / 2, ny / 2 + 16)
  File.open("grid.csv", "w") { |f| f.puts data.reverse }
end

main
puts `py plot_cfd.py -s`

__END__
