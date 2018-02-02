(defproject farg/acclivation "0.1.0-SNAPSHOT"
  :description "Acclivation laboratory"
  :url "https://github.com/bkovitz/acclivation"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[better-cond "1.0.1"]
                 [clj-time "0.14.2"]
                 [com.rpl/specter "1.1.0"]
                 [farg/util "0.1.0-SNAPSHOT"]
                 [farg/pmatch "0.1.0-SNAPSHOT"]
                 [farg/with-state "0.0.1-SNAPSHOT"]
                 [farg/x "0.1.0-SNAPSHOT"]
                 [incanter "1.5.7"]
                 [org.clojure/clojure "1.9.0"]
                 [org.clojure/math.numeric-tower "0.0.4"]
                 [org.clojure/math.combinatorics "0.1.1"]
                 [org.clojure/tools.macro "0.1.2"]
                 [org.clojure/tools.namespace "0.2.11"]
                 [org.clojure/tools.trace "0.7.9"]
                 [org.clojure/core.async "0.3.443"]
                 [net.mikera/core.matrix "0.61.0"]
                 [popen "0.3.1"]
                 [potemkin "0.4.4"]
                 [seesaw "1.4.5"]
                 [ubergraph "0.4.0"]]
  :main ^:skip-aot farg.acclivation
  :target-path "target/%s"
  :profiles {:uberjar {:aot :all}})
