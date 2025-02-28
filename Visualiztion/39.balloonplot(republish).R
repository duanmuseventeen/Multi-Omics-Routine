library(gplots)
library(datasets)
data(Titanic)

dframe <- as.data.frame(Titanic) # convert to 1 entry per row format
attach(dframe)
balloonplot(x=Class, y=list(Survived, Age, Sex), z=Freq, sort=TRUE)

# colorize: surviors lightblue, non-survivors: grey
Colors <- Titanic
Colors[,,,"Yes"] <- "skyblue"
Colors[,,,"No"] <- "grey"
colors <- as.character(as.data.frame(Colors)$Freq)

balloonplot(x=list(Age,Sex),
            y=list(Class=Class,
                   Survived=reorder.factor(Survived,new.order=c(2,1))
            ),
            z=Freq,
            zlab="Number of Passengers",
            sort=TRUE,
            dotcol = colors,
            show.zeros=TRUE,
            show.margins=TRUE)