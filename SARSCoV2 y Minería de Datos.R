##############
# Inicialización

install.packages("seqinr")
library(seqinr)
setwd("C:/Users/Elaia/Desktop/EVMHAD/")

secuencias <- read.fasta("Secuencias.fasta",forceDNAtolower = FALSE,
                         set.attributes = FALSE)
metaSecuencias <- read.csv("Secuencias.csv")


##############
# Valores iniciales

tamanioSecuencia   <- 1273
cantidadSecuencias <- length(secuencias)
#cantidadSecuencias <- 1000
listaDeAminoacidos <- c("A","R","N","D","C","Q","E","G","H","I",
                        "L","K","M","F","P","S","T","W","Y","V")
numeroAminoacidos <- length(listaDeAminoacidos)


##############
# Limpieza de datos

toleranciaAmbiguedad <- 0.02
maximoAmbiguedad <- floor(toleranciaAmbiguedad*tamanioSecuencia)
i <- 1
while(i<=cantidadSecuencias){
  aux <- 0
  for(j in 1:tamanioSecuencia){
    if(secuencias[[i]][j]=='X'){
      aux <- aux+1
    }
  }
  if(aux>maximoAmbiguedad){
    # Aqui me quita el dato "ambiguo"
    metaSecuencias <- metaSecuencias[metaSecuencias$Accession!=metaSecuencias$Accession[i],]
    secuencias[i] <- NULL
    cantidadSecuencias <- cantidadSecuencias-1
  }else{
    i <- i+1
  }
  rm(aux)
}
rm(i,j)
#Guardar proceso
#write.fasta(secuencias, names=names(secuencias), file.out="secuenciasPostLimpieza.fasta")
#write.csv(metaSecuencias, file="SecuenciasPostLimpieza.csv")


##############
# Secuencia consenso

matrizAux <- matrix(0,nrow=tamanioSecuencia,ncol=numeroAminoacidos)
colnames(matrizAux)<- listaDeAminoacidos
for(i in 1:cantidadSecuencias){
  for(j in 1:tamanioSecuencia){
    for(k in 1:numeroAminoacidos){
      if(secuencias[[i]][j]==listaDeAminoacidos[k]){
        matrizAux[j,k] <- matrizAux[j,k]+1
        break
      }
    }
  }
}
rm(i,j,k)
sum(rowSums(matrizAux)==cantidadSecuencias)
for(i in 1:tamanioSecuencia){       # Para poner las X en el consenso manualmente
  aux <- cantidadSecuencias-sum(matrizAux[i,])
  index <- which(matrizAux[i,]==max(matrizAux[i,]))
  matrizAux[i,index] <- aux + matrizAux[i,index]
}
sum(rowSums(matrizAux)==cantidadSecuencias)
rm(aux,index,i)
cadenaConsenso <- c()
for (i in 1:tamanioSecuencia) {
  for (j in 1:numeroAminoacidos) {     # Se podria resumir con un which()
    if (matrizAux[i,j]==max(matrizAux[i,])){
      cadenaConsenso <- c(cadenaConsenso,listaDeAminoacidos[j])
    }
  }
}
rm(i,j)
#Guardar proceso
#write.csv(matrizAux,"matrizAux.csv")


##############
# Revisar las mutaciones por posición

mutacionesPorPosicion <- c()
for (i in 1:tamanioSecuencia) {
  mutacionesPorPosicion <- c(mutacionesPorPosicion, sum(matrizAux[i,])-max(matrizAux[i,]))
}
#Guardar proceso
#write(mutacionesPorPosicion, "mutaciones.txt")
par(mfrow=c(1,1))
plot(mutacionesPorPosicion,main="USA")
logMutacionesPorPosicion <- log(mutacionesPorPosicion)
plot(logMutacionesPorPosicion,main="USA")



##############
# Filtrar las mutaciones por estado

# Estado 1
estado1 <- 'MI'
filtro1 <- which(metaSecuencias$USA[1:cantidadSecuencias]==estado1)
mutacionesPorPosicionEstado1 <- rep(0,tamanioSecuencia)
for(i in filtro1){
  for(j in 1:tamanioSecuencia){
    if(secuencias[[i]][j]!=cadenaConsenso[j]){
      mutacionesPorPosicionEstado1[j] <- mutacionesPorPosicionEstado1[j]+1
    }
  }
}
logMutacionesPorPosicionEstado1 <- log(mutacionesPorPosicionEstado1)

# Estado 2
estado2 <- 'CA'
filtro2 <- which(metaSecuencias$USA[1:cantidadSecuencias]==estado2)
mutacionesPorPosicionEstado2 <- rep(0,tamanioSecuencia)
for(i in filtro2){
  for(j in 1:tamanioSecuencia){
    if(secuencias[[i]][j]!=cadenaConsenso[j]){
      mutacionesPorPosicionEstado2[j] <- mutacionesPorPosicionEstado2[j]+1
    }
  }
}
logMutacionesPorPosicionEstado2 <- log(mutacionesPorPosicionEstado2)


par(mfrow=c(1,2))
plot(logMutacionesPorPosicionEstado1,main=estado1)
plot(logMutacionesPorPosicionEstado2,main=estado2)

par(mfrow=c(1,1))
plot(logMutacionesPorPosicionEstado1, col="red",main=estado1,
     ylim=c(0,log(max(length(filtro1),length(filtro2)))))
lines(logMutacionesPorPosicionEstado2, col="blue",type="p",main=estado2,
      ylim=c(0,log(max(length(filtro1),length(filtro2)))))

length(filtro1)
length(filtro2)


##############
# Normalización

muestrasEstado1 <- length(filtro1)
muestrasEstado2 <- length(filtro2)
variabilidadEstado1 <- logMutacionesPorPosicionEstado1/log(muestrasEstado1)
variabilidadEstado2 <- logMutacionesPorPosicionEstado2/log(muestrasEstado2)

par(mfrow=c(1,2))
plot(variabilidadEstado1,ylim=c(0,1),main=estado1)
plot(variabilidadEstado2,ylim=c(0,1),main=estado2)
