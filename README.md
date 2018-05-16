The codes in this repo were used to generate the results for our paper (add arxiv link of the paper)



### Folder Description

```
./ann_arbor : generate {truck+drone} path for the city of Ann Arbor, MI, USA
./gen_instances : generate test instances {demands, warehouse locations} for large scale simulation over Michigan, USA
./michigan : enerate {truck+drone} path for the southern part of the Michigan State, USA
./lib : shared scripts
```
See run.m in each folder to learn more detail on its contents

### Authors

* **Hyongju Park** 
* **Jinsun Liu**


### License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


### Acknowledgement

* "plot_google_map.m" by Zohar Bar-Yehuda (https://github.com/zoharby/plot_google_map) was used in some part of the code to display google map images
* "lldistkm.m" by M. Sohrabinia (https://www.mathworks.com/matlabcentral/fileexchange/38812-latlon-distance?focused=5250973&tab=function) was used in some part of the code to convert the shortest distance between any two GPS locations to kilometer
* "exchange2.m" by Jonas Lundgren (https://www.mathworks.com/matlabcentral/fileexchange/35178-tspsearch) was used to improve the solution of approximate TSP by 2-exchange heuristics
