##############
# Inicialización

install.packages("seqinr")
library(seqinr)
setwd("C:/Users/Elaia/Desktop/EVMHAD/")

secuencias <- read.fasta("Secuencias.fasta",as.string = TRUE,
                         forceDNAtolower = FALSE,set.attributes = FALSE)
metaSecuencias <- read.csv("Secuencias.csv")


##############
# Valores iniciales

tamanioSecuencias <- 1273
#cantidadSecuencias <- length(secuencias)
cantidadSecuencias <- 1000
listaDeAminoacidos <- c("A","R","N","D","C","Q","E","G","H","I",
                        "L","K","M","F","P","S","T","W","Y","V")
numeroAminoacidos <- length(listaDeAminoacidos)


##############
# Limpieza de datos

toleranciaAmbiguedad <- 0.02
maximoBasura <- floor(toleranciaAmbiguedad*tamanioSecuencias)
i <- 1
while(i<=cantidadSecuencias){
  aux <- 0
  for(j in 1:tamanioSecuencias){
    if(substring(secuencias[[i]],j,j)=='X'){
      aux <- aux+1
    }
  }
  if(aux>maximoBasura){
    # Aqui me quita el dato "basura"
    metaSecuencias <- metaSecuencias[metaSecuencias$Accession!=metaSecuencias$Accession[i],]
    secuencias[i] <- NULL
    cantidadSecuencias <- cantidadSecuencias-1
  }else{
    i <- i+1
  }
  rm(aux)
}
rm(i,j)


##############
# Secuencia consenso

matrizAux <- matrix(0,nrow=tamanioSecuencias,ncol=numeroAminoacidos)
colnames(matrizAux)<- listaDeAminoacidos
for(i in 1:cantidadSecuencias){
  for(j in 1:tamanioSecuencias){
    for(k in 1:numeroAminoacidos){
      if(substring(secuencias[[i]],j,j)==listaDeAminoacidos[k]){
        matrizAux[j,k] <- matrizAux[j,k]+1
        break
      }
    }
  }
}
sum(rowSums(matrizAux)==cantidadSecuencias)
for(i in 1:tamanioSecuencias){       # Para poner las X en el consenso manualmente
  aux <- cantidadSecuencias-sum(matrizAux[i,])
  index <- which(matrizAux[i,]==max(matrizAux[i,]))
  matrizAux[i,index] <- aux + matrizAux[i,index]
}
sum(rowSums(matrizAux)==cantidadSecuencias)
rm(aux, index)
cadenaConsenso <- c()
for (i in 1:tamanioSecuencias) {
  for (j in 1:numeroAminoacidos) {
    if (matrizAux[i,j]==max(matrizAux[i,])){
      cadenaConsenso <- c(cadenaConsenso,listaDeAminoacidos[j])
    }
  }
}


##############
# Revisar las mutaciones por posición

mutacionesPorPosicion <- c()
for (i in 1:tamanioSecuencias) {
  mutacionesPorPosicion <- c(mutacionesPorPosicion, sum(matrizAux[i,])-max(matrizAux[i,]))
}
plot(mutacionesPorPosicion)
logMutacionesPorPosicion <- log(mutacionesPorPosicion)
plot(logMutacionesPorPosicion)


##############
# Filtrar las mutaciones por estado

estado1 <- 'CA'
filtro1 <- which(metaSecuencias$USA[1:cantidadSecuencias]==estado1)
matrizAuxEstado1 <- matrix(0,nrow=tamanioSecuencias,ncol=numeroAminoacidos)
colnames(matrizAuxEstado1)<- listaDeAminoacidos
for(i in filtro1){
  for(j in 1:tamanioSecuencias){
    for(k in 1:numeroAminoacidos){
      if(substring(secuencias[[i]],j,j)==listaDeAminoacidos[k]){
        matrizAuxEstado1[j,k] <- matrizAuxEstado1[j,k]+1
        break
      }
    }
  }
}
mutacionesPorPosicionEstado1 <- c()
for (i in 1:tamanioSecuencias) {
  mutacionesPorPosicionEstado1 <- c(mutacionesPorPosicionEstado1, sum(matrizAuxEstado1[i,])-max(matrizAuxEstado1[i,]))
}
logMutacionesPorPosicionEstado1 <- log(mutacionesPorPosicionEstado1)
estado2 <- 'WA'
filtro2 <- which(metaSecuencias$USA[1:cantidadSecuencias]==estado2)
matrizAuxEstado2 <- matrix(0,nrow=tamanioSecuencias,ncol=numeroAminoacidos)
colnames(matrizAuxEstado2)<- listaDeAminoacidos
for(i in filtro2){
  for(j in 1:tamanioSecuencias){
    for(k in 1:numeroAminoacidos){
      if(substring(secuencias[[i]],j,j)==listaDeAminoacidos[k]){
        matrizAuxEstado2[j,k] <- matrizAuxEstado2[j,k]+1
        break
      }
    }
  }
}
mutacionesPorPosicionEstado2 <- c()
for (i in 1:tamanioSecuencias) {
  mutacionesPorPosicionEstado2 <- c(mutacionesPorPosicionEstado2, sum(matrizAuxEstado2[i,])-max(matrizAuxEstado2[i,]))
}
logMutacionesPorPosicionEstado2 <- log(mutacionesPorPosicionEstado2)

par(mfrow=c(1,2))
plot(logMutacionesPorPosicionEstado1)
plot(logMutacionesPorPosicionEstado2)

plot(logMutacionesPorPosicionEstado1, col="red",ylim=c(0,log(length(filtro2))))
lines(logMutacionesPorPosicionEstado2, col="blue", type = "p")

length(filtro1)
length(filtro2)


##############
# Normalización

muestrasEstado1 <- length(filtro1)
muestrasEstado2 <- length(filtro2)
variabilidadEstado1 <- logMutacionesPorPosicionEstado1/log(muestrasEstado1)
variabilidadEstado2 <- logMutacionesPorPosicionEstado2/log(muestrasEstado2)

par(mfrow=c(1,2))
plot(variabilidadEstado1,ylim=c(0,1))
plot(variabilidadEstado2,ylim=c(0,1))
