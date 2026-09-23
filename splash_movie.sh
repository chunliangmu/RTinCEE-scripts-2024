#!/usr/bin/bash

for i in no-apr--np21504 no-apr--np172032 merge2 merge4 merge12-com; do
for j in lumo fall; do
	echo && echo && echo ">>>>>> $i.$j <<<<<<" && echo && echo
	splash --movie $i/${j}_????0 -x 1 -y 2 -r 6
	rm -f splash.mp4
	ffmpeg -framerate 6 -i splash_%04d.png -r 60 -vb 50M -bt 100M -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -y $i.$j.splash.mp4
	tar -czvf $i.$j.splash.tgz splash*
	cp splash_0100.png $i.$j.splash_0100.png
	cp splash_0200.png $i.$j.splash_0200.png
	rm -f splash_????.png
done; done
