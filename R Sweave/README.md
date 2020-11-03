Die Auswertungen wurden mit R Sweave und der Satz mit Latex kombiniert durchgef�hrt. Die einzelnen Teile der Arbeit sind unter /parts, die Tabellen unter /tabellen und die Plots als PDF unter /plots zu finden. S�mtliche Auswertungen wurden in **Fehlende_Daten_pks.Rnw** in R durchgef�hrt (oder wurden mit \SweaveInput{...} daraus gestartet) und die Plots darin erzeugt. Die Tabellen passen sich bei �nderungen in der Auswertung automatisch an und sind so immer aktuell.

Die Daten, die durch die entsprechenden Simulationsabl�ufe ([...]_flow.R) erzeugt wurden sind in /data enthalten und werden von dort in **Fehlende_Daten_pks.Rnw** eingelesen.

F�r das Kompilieren sind unter Umst�nden zus�tzliche Schriftarten notwendig ("CM Roman") sowie Ghostscipt um diese in die PDF Plots und sp�ter in das durch PDFLatex erzeugte PDF einzubinden. (Es kann sonst auf anderen PCs zur Anzeige falscher Schriftarten kommen, wenn diese nicht dort lokal installiert sind.)
