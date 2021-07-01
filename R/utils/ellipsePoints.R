# helper function for plotting GPS observation errors
ellipsePoints = function(center, major, minor, angle, segments) {
  angles = (0:segments) * 2 * pi/segments
  unit.circle = cbind(cos(angles), sin(angles))
  cos.angle = cos(angle)
  sin.angle = sin(angle)
  Q = matrix(c(
    minor * cos.angle, minor * -sin.angle,
    major * sin.angle, major * cos.angle
  ), ncol = 2, byrow = TRUE)
  t(center + t(unit.circle %*% Q))
}