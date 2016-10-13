transform_to_0_1 = function(Xin) {
  if (is.numeric(Xin)) {(Xin - min(Xin)) / (max(Xin) - min(Xin))}
  else if (is.matrix(Xin)) {apply(Xin, 2, function(Xin2){(Xin2 - min(Xin2)) / (max(Xin2) - min(Xin2))})}
}
