#Dr. David Martínez Méndez
#Análisis de sistemas complejos discretos

library(BoolNet)
library(ggplot2)

#Carga de la red
net <- loadNetwork("tcell_activation_network.txt")

#Análisis general de la red
net
plotNetworkWiring(net)

#Cómputo y ploteo de atractores
attr <- getAttractors(net)
plotAttractors(attr)
getAttractorSequence(attr, 2)

#Cómputo y ploteo de atractores asincronos
attr2 <- getAttractors(net, type="asynchronous")

#Estados iniciales especificos
startState <- generateState(net, specs=c("TCR"=1, "CD28"=1, "CD8086"=1))
plotSequence(net, startState=startState)

#Cuencas de atracción
attr <- getAttractors(net, method="random", startStates=500)
plotStateGraph(attr, edge.arrow.size=0.01, piecewise=TRUE)

#Mutaciones
mutaciones <- fixGenes(net, c("TCR","CTLA4"), c(0,1))
attr <- getAttractors(mutaciones)
attr

#Análisis de robustez por ruido en transiciones
#Distancia de Hamming
r <- perturbTrajectories(net, measure="hamming", numSamples=100, flipBits=1)
r$value

#Distancia de atractores
r <- perturbTrajectories(net, measure="attractor", numSamples=100, flipBits=1)
r$value

#Análisis de robustez por ruido en estructura de la red
#Bitflip
perturbedNet <- perturbNetwork(net, perturb="functions", method="bitflip")
attr <- getAttractors(perturbedNet)
attr

#Shuffle
perturbedNet <- perturbNetwork(net, perturb="functions", method="shuffle")
attr <- getAttractors(perturbedNet)
attr
