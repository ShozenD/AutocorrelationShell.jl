function [RS,RD]=ACCorrCalc(R, w::OrthoFilter, L)


  P = Pfilter(w)
  Q = Qfilter(w)

  b = autocorr(P) / 2
  c = autocorr(Q) / 2

  b = vcat(b[end:-1:1], norm(P)^2, b)
  c = vcat(c[end:-1:1], norm(Q)^2, c)

  n = length(R)
  J = log2(n)

  RS = zeros(n, L+1)
  RD = zeros(n, L)
  RS[:,1] = R

  for j = 1:L
    RS[1:2^(J - j), j + 1] = subsample(ac_filter(RS[1:2^(J - j + 1), j]', b))';
    RD[1:2^(J - j), j] = subsample(ac_filter(RS[1:2^(J - j + 1), j]',c))';
  end
end
