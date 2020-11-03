Die Auswertungen wurden mit R Sweave und der Satz mit Latex kombiniert durchgeführt. Die einzelnen Teile der Arbeit sind unter /parts, die Tabellen unter /tabellen und die Plots als PDF unter /plots zu finden. Sämtliche Auswertungen wurden in **Fehlende_Daten_pks.Rnw** in R durchgeführt (oder wurden mit \SweaveInput{...} daraus gestartet) und die Plots darin erzeugt. Die Tabellen passen sich bei Änderungen in der Auswertung automatisch an und sind so immer aktuell.

Die Daten die durch die entsprechenden Simulationsabläufe ([...]_flow.R)erzeugt wurden sind in /data enthalten und werden von dort in **Fehlende_Daten_pks.Rnw** eingelesen.

Für das Kompilieren sind entsprechende Schriftarten notwendig ("CM Roman") sowie Ghostscipt um diese in die PDF Plots und später in das durch PDFLatex erzeugt PDF einzubinden. (Kann sonst auf anderen PCs zur Anzeige falscher Schriftarten führen, wenn diese nicht dort lokal installiert sind.)
