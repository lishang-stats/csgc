#' Calculate basic centred subgraph counts
#'
#' @param A adjacency matrix
#' @param P probability matrix
#' @param var.structure allow choose of bernoulli/poisson model.
#'
#' @return Unnormalised statistics (V.graph), variance of statistics (vr.graph),
#' normalised statistics (t.graph)
#' @export
#'
#' @examples
#' library(csgc)
#' k = 4
#' n = 100
#' prob = c(0.8, 0.1)
#' K = matrix(prob[2],k,k)
#' diag(K) = prob[1]
#' z = rep(1:k,each=n/k)
#' mat = gen_adj_sbm(K,z)
#' A = mat$A
#' P = mat$P
#' csgc(A,P,"bernoulli")
csgc = function(A, P, var.structure=c("bernoulli", "poisson")) {
  # INPUT
  # A     = adjacency matrix
  # P     = hypothesised connection probabilities
  # A and P need to have the same dimensions
  #
  # OUTPUT
  # A list with named vectors:
  # t:   normalised centred subgraph counts
  # V:   unnormalised centred subgraph counts
  # var: variances of V
  #
  # NOTE
  # csgc(A,0*A) will return the raw subgraph counts in V
  var.structure = match.arg(var.structure)

  if (var.structure == "bernoulli"){
    if (sum(A>1)/prod(dim(A)) > 0.01){
      warning("Results may be inaccurate: too many multiedges and selfedges.")
    }
    # Clean up A and P
    diag(A) = 0 # remove loops
    diag(P) = 0
    # ignore those entries where P>1
    w = which(P>1,arr.ind=TRUE)
    P[w] = 1
    # Calculate variance matrix
    vr = P*(1-P)
    # Calculate centred indicators
    Yh = A - P
    # Remove those where P>1
    Yh[w] = 0
  }

  if (var.structure == "poisson") {
    # Calculate variance matrix
    vr = P
    # Calculate centred indicators
    Yh = A - P
  }

  # Calculate some intermediary matrices
  Yh.2p = Yh %*% Yh
  Yh.3p = Yh.2p %*% Yh
  Yh.4p = Yh.3p %*% Yh
  Yh.5p = Yh.4p %*% Yh
  Yh.2 = Yh*Yh
  Yh.3 = Yh.2*Yh
  Yh.4 = Yh.3*Yh
  Yhp.2.1 = Yh.2 %*% Yh
  Yhp.2.2 = Yh.2 %*% Yh.2
  Yhp.2.2.1 = Yhp.2.2 %*% Yh
  Yhp.2.2p = Yh.2 %*% Yh.2p
  Yhp.2p.2 = Yh.2p %*% Yh.2
  Yh.2p.2p.1 = Yh.2p*Yh.2p*Yh
  Yhp.3.1 = Yh.3 %*% Yh
  Yhp.3.2p = Yh.3 %*% Yh.2p
  Yhr = rowSums(Yh)
  Yhr.2 = Yhr*Yhr
  Yhr.3 = Yhr.2*Yhr
  Yhr.4 = Yhr.3*Yhr
  Yh2r = rowSums(Yh.2)
  Yh2r.2 = Yh2r*Yh2r
  Yhrr = Yh2r*Yhr.2
  Yhd.3p = diag(Yh.3p)
  Yhd.3p.1 = Yhd.3p %*% Yh
  Yhd.3p.2 = Yhd.3p %*% Yh.2

  vr.2p = vr %*% vr
  vr.3p = vr.2p %*% vr
  vr.4p = vr.3p %*% vr
  vr.5p = vr.4p %*% vr
  vr.2 = vr*vr
  vr.3 = vr.2*vr
  vr.4 = vr.3*vr
  vrp.2.1 = vr.2 %*% vr
  vrp.2.2 = vr.2 %*% vr.2
  vrp.2.2.1 = vrp.2.2 %*% vr
  vrp.2.2p = vr.2 %*% vr.2p
  vrp.2p.2 = vr.2p %*% vr.2
  vr.2p.2p.1 = vr.2p*vr.2p*vr
  vrp.3.1 = vr.3 %*% vr
  vrp.3.2p = vr.3 %*% vr.2p
  vrr = rowSums(vr)
  vrr.2 = vrr*vrr
  vrr.3 = vrr.2*vrr
  vrr.4 = vrr.3*vrr
  vr2r = rowSums(vr.2)
  vr2r.2 = vr2r*vr2r
  vrrr = vr2r*vrr.2
  vrd.3p = diag(vr.3p)
  vrd.3p.1 = vrd.3p %*% vr
  vrd.3p.2 = vrd.3p %*% vr.2

  # 1
  V.twostar = sum((Yh.2p)[upper.tri(Yh)])
  vr.twostar = sum((vr.2p)[upper.tri(vr)])
  # 2
  V.triangle = trace(Yh.3p)/6
  vr.triangle = trace(vr.3p)/6
  # 3
  V.fourcycle = (trace(Yh.4p) - 4 * sum((Yhp.2.2)[upper.tri(Yh)]) - sum(Yh.4) )/8
  vr.fourcycle = (trace(vr.4p) - 4 * sum((vrp.2.2)[upper.tri(vr)]) - sum(vr.4) )/8
  # 4
  V.threepath = sum((Yh.3p)[upper.tri(Yh)]) + sum((Yh.3)[upper.tri(Yh)]) - sum(Yhp.2.1)
  vr.threepath = sum((vr.3p)[upper.tri(vr)]) + sum((vr.3)[upper.tri(vr)]) - sum(vrp.2.1)
  # 5
  V.threestar = (sum(Yhr.3) - 3*sum(Yhp.2.1) + 2*sum(Yh.3))/6
  vr.threestar = (sum(vrr.3) - 3*sum(vrp.2.1) + 2*sum(vr.3))/6
  # 6
  V.triangleapp = (sum(Yhd.3p.1) - 2*trace(Yhp.2p.2))/2
  vr.triangleapp = (sum(vrd.3p.1) - 2*trace(vrp.2p.2))/2
  # 7
  V.twotriangle = (sum(Yh.2p.2p.1) - trace(Yhp.2.2.1))/4
  vr.twotriangle = (sum(vr.2p.2p.1) - trace(vrp.2.2.1))/4
  # 8
  V.fivecycle = (trace(Yh.5p) - 5*sum(Yhd.3p.2) + 5*trace(Yhp.3.2p))/10
  vr.fivecycle = (trace(vr.5p) - 5*sum(vrd.3p.2) + 5*trace(vrp.3.2p))/10
  # 9
  V.fourpath = (2*sum((Yh.4p)[upper.tri(Yh)]) - 2*sum(Yhp.2.2p) +
                  2*sum(Yhp.3.1) + 3*sum(Yhp.2.2) + 3*trace(Yhp.2.2p) -
                  2*sum(Yh.4) - sum(Yhrr) - 2*sum(Yhd.3p.1))/2
  vr.fourpath = (2*sum((vr.4p)[upper.tri(vr)]) - 2*sum(vrp.2.2p) +
                   2*sum(vrp.3.1) + 3*sum(vrp.2.2) + 3*trace(vrp.2.2p) -
                   2*sum(vr.4) - sum(vrrr) - 2*sum(vrd.3p.1))/2
  # 10
  V.fourstar = (sum(Yhr.4) - 6*sum(Yhrr) + 3*sum(Yhp.2.2)
                + 8*sum(Yhp.3.1) - 6*sum(Yh.4))/24
  vr.fourstar = (sum(vrr.4) - 6*sum(vrrr) + 3*sum(vrp.2.2)
                 + 8*sum(vrp.3.1) - 6*sum(vr.4))/24

  V = c(V.twostar, V.triangle, V.fourcycle, V.threepath, V.threestar,
        V.triangleapp, V.twotriangle, V.fivecycle, V.fourpath, V.fourstar)
  names(V) = c("V.twostar", "V.triangle", "V.fourcycle", "V.threepath", "V.threestar",
               "V.triangleapp", "V.twotriangle", "V.fivecycle", "V.fourpath", "V.fourstar")
  vr = c(vr.twostar, vr.triangle, vr.fourcycle, vr.threepath, vr.threestar,
         vr.triangleapp, vr.twotriangle, vr.fivecycle, vr.fourpath, vr.fourstar)
  names(vr) = c("var.twostar", "var.triangle", "var.fourcycle", "var.threepath", "var.threestar",
                "var.triangleapp", "var.twotriangle", "var.fivecycle", "var.fourpath", "var.fourstar")
  t = V/sqrt(vr)
  names(t) = c("t.twostar", "t.triangle", "t.fourcycle", "t.threepath", "t.threestar",
               "t.triangleapp", "t.twotriangle", "t.fivecycle", "t.fourpath", "t.fourstar")
  ans = list(t=t, V=V, var=vr)
  return(ans)
}
