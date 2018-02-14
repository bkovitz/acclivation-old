(ns farg.acclivation
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.core.reducers :as r]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.math.combinatorics :as combo]
            [clojure.math.numeric-tower :as math]
            [clojure.java.io :as io :refer [file writer]]
            [clj-time.local :as ltime]
            [com.rpl.specter :as S]
            [farg.util :refer [dd dde rand choose choose-one with-rng-seed
                               mround defopts with-*out*]
             :as util]
            [farg.with-state :refer [with-state]]
            [farg.acclivation.sa :as sa]
            [farg.acclivation.hill :as hill]
            [potemkin :refer [fast-memoize] :rename {fast-memoize memoize}]
            [ubergraph.core :as uber])
  (:gen-class))

;;;; Naming conventions:
;;;;
;;;; w    A fitness function ("worth" in the population-biology literature).
;;;; ph   phenotype, a vector. Destructured into [x1 x2].


;NEXT
;
; fitness-plot-data   and plot the phenotype fitness  DONE
; genotype->phenotype  DONE
; mutate  DONE
; nextgen  DONE

; save pop in atom  DONE
; fitness function with a "ridge": changing optimal value for x1  DONE
; run-epoch  DONE
; vary optimal spot on ridge each epoch  DONE

(declare add-fitness)

#_(def ^:dynamic *pop-file* *out*) ;output file for population
(def ^:dynamic *data-dir* (io/file "."))  ;TODO rm

(defn matlab [command]
  (clojure.java.shell/sh "matlab" "-nodesktop" "-nosplash" "-r"
                         (str command "; exit;") ">>" "matlab.out" "2>&1"))

(def lastpop (atom nil))

(defn printpop
 ([] (printpop @lastpop))
 ([p]
  (run! println (:individuals p))
  (run! println (:history p))))

(defn default-to [x default]
  (if (some? x) x default))

(defn strip-type [x]
  (vary-meta x dissoc :type))

(defn clear-data-dir [d]
  (when (and d (not= d (io/file ".")))
    (util/rm-recursively d)))

(defn make-file [basename {:keys [data-dir] :as opts}]
  (let [file (io/file data-dir basename)]
    (io/make-parents file)
    file))

(defn task-k->file
  [{:keys [epoch generation]} suffix opts]
  (make-file (cond
               (and (nil? epoch) (nil? generation))
                 suffix
               (nil? generation)
                 (str "epoch" epoch suffix)
               (str "epoch" epoch "-gen" generation suffix))
             opts))

;;; edn ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn strip-type [x]
  (vary-meta x #(dissoc % :type)))

(defn edge->vec [g e]
  [(uber/src e) (uber/dest e) (uber/attr g e :weight)])

(defn graph->map [g]
  {:nodes (->> g uber/nodes vec)
   :edges (->> g uber/edges (map #(edge->vec g %)))})

(defn map->graph [m]
  (with-state [g (uber/digraph)]
    (doseq [node (:nodes m)]
      (uber/add-nodes node))
    (uber/add-edges* (:edges m))))

(defn genotype->edn [gt]
  (str "#farg.acclivation/genotype"
       (-> gt
           strip-type
           (update :graph graph->map)
           (dissoc :w)
           prn-str)))

(defn map->genotype [m]
  (-> m
      (update :graph map->graph)
      (with-meta {:type :acclivation})))

(def edn-readers {'farg.acclivation/genotype map->genotype})

(defn edn-read [pushback-reader]
  (clojure.edn/read {:readers edn-readers, :eof :eof} pushback-reader))

(defn edn-slurp [file]
  (let [pbr (java.io.PushbackReader. (io/reader file))]
    (->> (repeatedly #(edn-read pbr))
      (take-while #(not= :eof %)))))
;  (clojure.edn/read {:readers edn-readers, :eof :eof}
;    (java.io.PushbackReader. (io/reader file))))

;;; phenotype fitness ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-inverted-v [target radius]
  (fn [x]
    (let [dist (Math/abs (- x target))]
      (if (> dist radius)
          0.0
          (let [proportion (/ (- radius dist)
                              radius)]
            (* proportion proportion))))))

(defn w-equal [[x1 x2]]
  (let [invv (make-inverted-v 0.0 0.2)]
    (invv (Math/abs (- x1 x2)))))

(def w-ratio
  (let [invv (make-inverted-v 1.50 0.1)]
    (fn [[x1 x2]]
      (if (zero? x1)
          0.0
          (invv (/ x2 x1))))))

(def w2
  (let [near1 (util/piecewise-linear 0.0 0.0, 1.0 1.0, 1.0 0.0)]
    (fn [[x1 x2]]
      (+ (near1 x1) (near1 x2)))))

(defn w12 [ph]
  (* (w-ratio ph) (w2 ph)))

(def some-points
  (with-rng-seed 1
    (vec (repeatedly 20 #(vector (rand -1.0 +1.0) (rand -1.0 +1.0))))))

(defn w-many-small-hills [[x1 x2]]
  (* (Math/cos (* 20 x1))
     (Math/sin (* 20 x2))))

(defn w-distance [[x1 x2]]
  (let [center [-0.43 0.67]
        invv (make-inverted-v 0.0 0.5)]
    #_(* 5.0 (invv (util/distance [x1 x2] center)))
    (- 1.0 (util/distance [x1 x2] center))
    ))

#_(defn w [ph]
  (+ (w12 ph) (w-many-small-hills ph) (w-distance ph)))

#_(defn w [ph]
  (if (some #(< (Math/abs %) 0.02) ph)
    0.0
    (+ (w-many-small-hills ph)
       (* 9.0 (w12 ph) (w-distance ph)))))

#_(defn w [ph]
  (+ (w-many-small-hills ph) (w-distance ph)))

#_(defn w [ph]
  (* #_(w-equal ph) (Math/pow (w-distance ph) 2.0)))

(def w-equal
  (let [invv (make-inverted-v 0.0 0.5)]
    (fn [[x1 x2]]
      (invv (- x1 x2)))))

(defn make-w [c]
  (fn [ph]
    (+ (w-many-small-hills ph)
       (- 1.0 (util/distance [c c] ph)))))

(defn make-w [c]
  (fn [[x1 x2 :as ph]]
    (+ (w-many-small-hills ph)
       (* 5.0 (w-equal ph) (- 1 (Math/abs (- x1 c)))))))

(defn make-epoch-w
  "Returns [w w-source] where w is a randomly generated phenotype fitness
  function and w-source is the source code that generates it."
  []
  (let [c (rand -1.0 +1.0)
        w-source `(make-w ~c)
        w (eval w-source)]
    [w w-source]))

(defn wfn [{:keys [w w-source]}]
  (assert (some? w-source))
  (if (some? w) w (memoize (eval w-source))))

;;; genotype fitness ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn extract-phenotype [g]
  [(sa/a g :p1) (sa/a g :p2)])

(defn spread-activation [g]
  (-> (:graph g)
      (sa/spread-activation (zipmap [:g1 :g2] (:numbers g)) :iterations 20)))

(defn genotype->phenotype [g]
  (-> (spread-activation g)
      (extract-phenotype)))

;(def genotype->phenotype (memoize genotype->phenotype))

(def ^:dynamic *genotype->phenotype* genotype->phenotype)

(defn gfn
  "Genotype gt's genotype->phenotype function, as a function of :numbers."
  [gt]
  (fn [xx]
    (let [gt' (assoc gt :numbers xx)]
      (*genotype->phenotype* gt'))))

(defn vfn
  "Return's genotype gt's virtual fitness function."
  [gt]
  (let [xx->ph (gfn gt)
        w (wfn gt)]
    (fn [xx]
      (-> xx xx->ph w))))

(defn penalize-zeros [w ph]
  (if (some zero? ph)
    -10.0
    (w ph)))

(defn genotype-fitness [g]
  (if-let [fitness (:fitness g)]
    fitness
    (let [w (wfn g)
          _ (assert (some? w) "Need phenotype fitness function")
          w #(penalize-zeros w %)]
      (->> (*genotype->phenotype* g)
           (w)))))

(defn add-phenotype [g]
  (if (contains? g :phenotype)
      g
      (assoc g :phenotype (*genotype->phenotype* g))))

(defn add-fitness [g]
  (let [g (add-phenotype g)]
    (if (contains? g :fitness)
        g
        (assoc g :fitness (genotype-fitness g)))))

;;; making genotypes ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-node [n]
  {:pre [(integer? n)]}
  (keyword (str \n n)))

(defn random-weight []
  (choose-one -1.0 1.0))

(defn make-random-graph []
  (let [n-extra-nodes (util/rand-int 1 5)
        extra-nodes (->> (take n-extra-nodes (iterate inc 1))
                         (map make-node))
        nodes (into [:g1 :g2 :p1 :p2] extra-nodes)
        n-edges (util/rand-int 5 (* 5 (count nodes)))]
    (with-state [g (apply uber/digraph nodes)]
      (dotimes [n n-edges]
        (bind [from to] (util/choose-with-replacement 2 nodes))
        (bind weight (random-weight))
        (uber/add-edges [from to weight]))
      (uber/add-edges [:g1 (choose extra-nodes) (random-weight)])
      (uber/add-edges [:g2 (choose extra-nodes) (random-weight)])
      (uber/add-edges [(choose extra-nodes) :p1 (random-weight)])
      (uber/add-edges [(choose extra-nodes) :p2 (random-weight)])
      )))

(def empty-genotype
  ^{:type :acclivation}
  {:type :genotype})

(defn make-random-genotype []
  (assoc empty-genotype
     :numbers [(util/rand -1.0 1.0) (util/rand -1.0 1.0)]
     :graph (make-random-graph)))

(def empty-population
  ^{:type :acclivation}
  {:type :population
   :generation 0
   :individuals nil
   :history []
   :epoch 1})

(defn add-basic-data [population gt]
  (assert (:w population))
  (-> gt
    (merge (select-keys population [:w-source :w :epoch :generation]))
    (add-phenotype)
    (add-fitness)))

(defn add-basic-data-to-everybody [population]
  (update population :individuals (fn [individuals]
                                    (map #(add-basic-data population %)
                                         individuals))))

(defn make-random-population [population n]
  (assoc population
         :individuals (->> (repeatedly make-random-genotype)
                           (map #(add-basic-data population %))
                           (take n)
                           vec)
         :population-size n))

;;; mutation ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn clamp [x]
  (util/clamp [-1.0 1.0] x))

(defn turn-knob [g]
  (let [i (choose [0 1])
        Δ (util/rand -0.02 0.02)]
    (update-in g [:numbers i] #(sa/squash (+ % Δ)))))

(defn turn-knob [g]
  (let [i (choose [0 1])
        Δ (* 0.01 (util/sample-normal))]
    (update-in g [:numbers i] #(clamp (+ % Δ)))))

(defn n-node? [node]
  {:pre [(keyword? node)]}
  (= \n (first (name node))))

(defn has-n-node? [g n]
  {:pre [(integer? n)]}
  (uber/has-node? (:graph g) (make-node n)))

(defn n-nodes-of [g]
  (->> (uber/nodes (:graph g))
       (remove (comp not n-node?))))

(defn first-available-node [g]
  (let [n (first (drop-while #(has-n-node? g %) (iterate inc 1)))]
    (make-node n)))

(defn add-node [g]
  (let [new-node (first-available-node g)
        graph (with-state [graph (:graph g)]
                (bind nodes (uber/nodes graph))
                (uber/add-nodes new-node)
                (uber/add-edges [new-node (choose nodes) (random-weight)])
                (uber/add-edges [(choose nodes) new-node (random-weight)]))]
    (assoc g :graph graph)))

(defn remove-node [g]
  (update g :graph #(uber/remove-nodes % (choose (n-nodes-of g)))))

(defn add-edge [g]
  (update g :graph #(let [nodes (uber/nodes %)
                          from (choose nodes)
                          to (choose nodes)]
                      (uber/add-edges % [from to (random-weight)]))))

(defn remove-edge [g]
  (let [graph (:graph g)
        edge (choose (uber/edges graph))]
    (if edge
      (assoc g :graph (uber/remove-edges graph edge))
      g)))
; (update g :graph #(uber/remove-edges % (choose (uber/edges %)))))

(defn move-edge [g]
  (let [graph (:graph g)
        nodes (uber/nodes graph)]
    (if-let [e (choose (uber/edges graph))]
      (let [weight (uber/weight graph e)]
        (assoc g :graph (with-state [graph graph]
          (bind from (choose nodes))
          (bind to (choose nodes))
          (uber/remove-edges e)
          (uber/add-edges [from to weight]))))
      (add-edge g))))

(def mutations
  [[turn-knob 10]
   [move-edge 1]
   [add-node 1]
   [remove-node 1]
   [add-edge 1]
   [remove-edge 1]])

(defn mutate [g]
  ((util/weighted-choice mutations) g))

(defn mutate-n [n g]
  (take n (repeatedly #(mutate g))))

(defn mutate-population [n p]
  (assoc p :individuals (->> (:individuals p)
                             (mapcat #(mutate-n n %)))))

;;; crossover ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn take-first-half [coll]
  (let [coll (seq coll)]
    (take (/ (count coll) 2) coll)))

(defn take-second-half [coll]
  (let [coll (seq coll)]
    (drop (/ (count coll) 2) coll)))

#_(defn crossover [g1 g2]
  (let [graph1 (:graph g1), graph2 (:graph g2)
        graph (-> (uber/digraph)
                  (uber/add-nodes* (uber/nodes graph1))
                  (uber/add-nodes* (uber/nodes graph2))
                  (uber/add-edges* (uber/edges graph1))
                  (uber/add-edges* (uber/edges graph2)))
        n1 ((:numbers g1) 0)
        n2 ((:numbers g2) 1)]
    [(assoc empty-genotype
            :numbers [n1 n2]
            :graph graph)]))

(defn about-half-of [coll]
  (keep #(when (< (rand) 0.7) %) coll))

;(def about-half-of take-first-half)

(defn mandatory-nodes [g]
  (->> (uber/nodes g)
       (filter #(#{\g \p} (-> % name str first)))))

#_(defn crossover [g1 g2]
  (let [graph1 (:graph g1), graph2 (:graph g2)
        graph (-> (uber/digraph)
                  (uber/add-nodes* (uber/nodes graph1))
                  (uber/add-nodes* (uber/nodes graph2))
                  (uber/add-edges* (uber/edges graph1))
                  (uber/add-edges* (uber/edges graph2)))]
    [(assoc empty-genotype
            :numbers (:numbers (util/choose-one g1 g2))
            :graph graph)]))

(defn crossover [g1 g2]
  (let [graph1 (:graph g1), graph2 (:graph g2)
        graph (-> (uber/digraph)
                  (uber/add-nodes* (mandatory-nodes graph1))
                  (uber/add-nodes* (uber/nodes graph1))
                  (uber/add-nodes* (uber/nodes graph2))
                  (uber/add-edges* (about-half-of (uber/edges graph1)))
                  (uber/add-edges* (about-half-of (uber/edges graph2))))]
    [(assoc empty-genotype
            :numbers (:numbers g1)
            :graph graph)]))

#_(defn crossover [g1 g2]
  [(assoc empty-genotype
          :numbers (:numbers g1)
          :graph (:graph g2))])

(defn tournament-selection [tourney-size fitness-func individuals]
  (apply max-key fitness-func
    (util/choose-without-replacement tourney-size individuals)))

(defn random-pair [coll]
  (util/choose-without-replacement 2 coll))

(defn mutate-and-crossover
  [tourney-size mutate crossover n-by-mutation n-by-crossover p]
  (let [parents (:individuals p)
        biased-choice
          #(tournament-selection tourney-size genotype-fitness parents)
        mutants (->> (repeatedly biased-choice)
                     (map mutate)
                     (mapcat concat)
                     distinct
                     (take n-by-mutation)
                     vec)
        crosses (->> (repeatedly #(vector (biased-choice) (biased-choice)))
                     (map #(apply crossover %))
                     (mapcat concat)
                     distinct
                     (take n-by-crossover))]
    (assoc p :individuals (into mutants crosses))))

#_(defn tournament-selection-on-population [f-fitness n p]
  (loop [winners #{}, n-tries 0]
    (if (or (>= (count winners) (:population-size p))
            (>= n-tries (* 2 (:population-size p))))
        (assoc p :individuals winners)
        (recur (conj winners (tournament-selection f-fitness n p))
               (inc n-tries)))))

(defn next-generation [vary select population]
  (-> population
      vary
      add-basic-data-to-everybody
      select
      #_(update :generation inc)))

;;; dot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn dot-edge [graph edge]
  (str (name (:src edge)) " -> " (name (:dest edge)) \space
       "[label=" (mround (uber/weight graph edge)) "]; "))
       ;"[label=" (mround (sa/a graph edge)) "]; "))

(defn dot-node [graph node]
  (str (name node) \space "[label=" (mround (sa/a graph node)) "] "))

(defn make-dot [graph]
  (str "digraph g {
{ rank=source edge [style=\"invis\"] g1 -> g2 }
{ rank=sink edge [style=\"invis\"] p1 -> p2 }
"
       (apply str (map #(dot-node graph %) (uber/nodes graph)))
       (apply str (map #(dot-edge graph %) (uber/edges graph)))
       "}"))

(defn dot [genotype]
  (->> genotype spread-activation make-dot))

(defn print-genotype-graph [genotype]
  (println (dot genotype)))

;;; print-method ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn graph-cstr [graph]
  (str (sort (uber/nodes graph)) \space
       (uber/count-edges graph) " edges"))

(defn genotype-cstr [g]
  (str (-> (strip-type g)
           (select-keys [:numbers :phenotype :fitness])) \space
       (graph-cstr (:graph g))))

(defn population-cstr [p]
  (str (strip-type (select-keys p [:generation :epoch])) \newline
       (count (:individuals p))))

(defn convenient-str ^String [x]
  (if (map? x)
      (case (:type x)
        :genotype (genotype->edn x) ;(genotype-cstr x)
        :population (population-cstr x)
        (apply str (strip-type x)))
      (apply str (strip-type x))))

(defmethod print-method :acclivation [v ^java.io.Writer w]
  (.write w (convenient-str v)))

(defn genotype-str [g]
  (->> (add-fitness g)
       (genotype-cstr)))

;;; printing ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

#_(defn print-population [p]
  (let [p (update p :individuals #(map add-fitness %))
        individuals (sort-by :fitness (:individuals p))]
    (with-*out* *pop-file*
      (when (not (zero? (:generation p)))
        (println))
      (println (strip-type (select-keys p [:generation :epoch])))
      (println (clojure.string/join \newline (map genotype-str individuals)))
      (print-genotype-graph (last individuals)))))

(defn best-of [p]
  (let [p (update p :individuals #(map add-fitness %))
        individuals (sort-by :fitness (:individuals p))]
    (last individuals)))
;    (->> (last individuals)
;         (add-basic-data p))))

(defn best-fitness-of [p]
  (:fitness (best-of p)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn trim-round-off-error [x]
  (/ (Math/round (* 100000000.0 x))
     100000000.0))

(defn normalize [x]
  (if (zero? x)
    0.0  ;force 0.0; prevent -0.0
    (-> x (* 100000) (Math/round) (/ 100000.0))))

(def normalized-range
  (->> (range -1.0 +1.0000000001 0.01)
       (map normalize)
       vec))

(defn print-fn [f]
  (->> (for [x normalized-range, y normalized-range]
         [x y])
       (pmap (fn [[x y]]
               [x y (f [x y])]))
       (run! #(apply println %))))

(defn save-fn [f filename]
  (with-*out* (io/writer filename)
    (print-fn f)))

(defn print-range
 ([gt]
  (print-range (gfn gt) (wfn gt)))
 ([gfn wfn]
  (->> (for [x normalized-range, y normalized-range]
         [x y])
       (pmap (fn [[x y]]
               (let [[phx phy] (gfn [x y])]
                 [phx phy (wfn [phx phy]) x y])))
       (run! #(apply println %)))))

(defn save-range
 ([gt filename]
  (save-range (gfn gt) (wfn gt) filename))
 ([gfn wfn filename]
  (with-*out* (io/writer filename)
    (print-range gfn wfn))))

;IDEA What's the average fitness?
;IDEA Look at where the fitness goes varying only one axis at a time.
;See the fitness function of x, holding y constant.
;IDEA See if doing a few mutations at once works better.
;IDEA Make a fitness landscape with concentric moats.
;IDEA See if a genome a single-dimensional start vector can find the diagonal
;ridge.

(defn make-vdata []
  (let [invv (make-inverted-v 5.0 2.0)]
    (doseq [x (range 0.0 10.0 0.05)]
      (println x (invv x)))))

#_(defn print-w
  "Prints the fitness function so we can compare it to the virtual fitness
  functions in the genotypes."
  [& {:keys [step] :or {step 0.005}}]
  (let [xrange (range -1.0 1.000000000001 step)
        yrange (range -1.0 1.000000000001 step)]
    (doseq [x xrange, y yrange]
      (let [x (trim-round-off-error x)
            y (trim-round-off-error y)
            fitness (w [x y])]
        (println x y fitness x y)))))

;;; measuring acclivity ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn acclivity
 ([f]
  (acclivity f 0.01 2 10))
 ([f step dimension n-climbers]
  (let [all-results (hill/run-climbers f step dimension n-climbers)]
    (run! println all-results)
    (util/average (map first all-results)))))

(defn save-vfn [genotype filename]
  (binding [*genotype->phenotype* genotype->phenotype]
    (save-fn (vfn genotype) filename)))

(defn save-gfn [genotype filename]
  (binding [*genotype->phenotype* genotype->phenotype]
    (save-fn (gfn genotype) filename)))

(defn genotype-acclivity
 ([genotype]
  (genotype-acclivity genotype 0.01 2 10))
 ([genotype step dimension n-climbers]
  (acclivity (vfn genotype) step dimension n-climbers)))

(defn fn-acclivity
  "Returns only average height reached by hill-climbers. Silent."
 ([f]
  (fn-acclivity f 0.01 2 10))
 ([f step dimension n-climbers]
  (->> (hill/run-climbers f step dimension n-climbers)
    (map first)
    util/average)))

(defn wfn-acclivity
 ([gt]
  (wfn-acclivity gt 0.01 2 10))
 ([gt step dimension n-climbers]
  (fn-acclivity (wfn gt) step dimension n-climbers)))

(defn vfn-acclivity
 ([gt]
  (vfn-acclivity gt 0.01 2 10))
 ([gt step dimension n-climbers]
  (fn-acclivity (vfn gt) step dimension n-climbers)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn now-suffix []
  (ltime/format-local-time (ltime/local-now) :basic-date-time-no-ms))

(defopts let-ga-opts
  {:keys [generations population-size n-mutants n-crosses f-mutate
                     f-cross tourney-size vary select fitness seed
                     radius step dimension n-epochs epochs-data]
   :or {generations 20, population-size 40, tourney-size 5,
        f-mutate (partial mutate-n 1), f-cross crossover,
        fitness genotype-fitness, seed 1, dimension 2, n-epochs 1,
        radius nil, step 0.01, epochs-data true}}
  n-mutants (default-to n-mutants (int (* population-size 0.71)))
  n-crosses (default-to n-crosses (- population-size n-mutants))
  vary (default-to vary
         (partial mutate-and-crossover
                  tourney-size f-mutate f-cross n-mutants n-crosses))
  ;select (partial tournament-selection-on-population fitness tourney-size)
  select identity
  next-gen (partial next-generation vary select)
  data-dir (io/file (str "seed" seed "-d" dimension))
;  best-genotype-file (io/file dirname "best-genotype")
;  best-fitness-file (io/file dirname "best-fitness")
;  hill-climbers-file (io/file dirname "hill-climbers")
;  generations-file (io/file dirname "generation")
;  prefix (str "seed" seed "-") ;prefix for output filenames
;  pop-file (file (str prefix "pop"))
;    ;population of each generation
;  vfit-file (file (str prefix "vfit"))
;    ;fitness-as-seen-by best of last gen
  )

#_(defn print-phenotype-fitness-fn [& opts]
  (let-ga-opts opts
    (with-*out* (writer "phenotype-fitness")
      (print-fitness-fn (fn [ph] [ph (w ph)])))))

#_(defn save-gen-data [{:keys [epoch generation] :as population}]
  (let [data-file (io/file *data-dir*
                          (str "best-genotype"
                               "-epoch" epoch
                               "-gen" generation))
        best-genotype (best-of population)]
    (io/make-parents data-file)
    (with-*out* (writer data-file)
      (prn best-genotype))
))

(defn save-gen-data [{:keys [epoch generation individuals] :as population}]
  (let [data-file (io/file *data-dir*
                          (str "epoch" epoch
                               "-gen" generation
                               ".pop"))]
    (io/make-parents data-file)
    (with-*out* (writer data-file)
      (run! prn individuals))))

(defn save-dot [dot-file gt opts]
  (with-open [dot-file (io/writer dot-file)]
    (with-*out* dot-file
      (print-genotype-graph gt))))

(defn prep-for-save [x]
  (if (= java.io.File (type x))
    (str x)
    x))

(defn save-opts [{:keys [data-dir] :as opts}]
  (let [opts-file (io/file data-dir "opts")
        opts (->> opts (keep (fn [[k v]]
                               (when (not (fn? v))
                                 [k (prep-for-save v)])))
                       (map vec)
                       (into {}))]
    (io/make-parents opts-file)
    (spit opts-file (pr-str opts))))

(defn accumulate-data [{:keys [generation] :as population}]
  (with-state [population population]
    add-basic-data-to-everybody
    (update :individuals #(map add-fitness %))
    (assoc :individuals (sort-by :fitness > (:individuals population)))
    ;(update :individuals #(map add-w-source population %))
    (update :history conj (best-fitness-of population))
    -- (reset! lastpop population)
    ))
  
(defn w-args [{:keys [w-source] :as gt}]
  (drop 1 w-source))

(defn n-edges [{:keys [graph] :as gt}]
  (-> graph uber/edges count))

(defn genotype-data [gt]
  [(:epoch gt) (:generation gt)
   (wfn-acclivity gt) (vfn-acclivity gt) (w-args gt)
   (:numbers gt) (n-edges gt) (:phenotype gt) (:fitness gt)])

(def genotype-data (memoize genotype-data))

(defn mround-floats [data]
  (S/transform [(S/walker float?)] util/mround data))

(defn mkfile [file]
  (if (= java.io.File (type file))
    file
    (io/file file)))

(defn readpop [filename]
  (-> filename edn-slurp vec))

(defn readbest
  "The best genotype is assumed to be the first genotype in the file."
  [filename]
  (-> filename edn-slurp first))

(defn parse-filename [file]
  (let [file (mkfile file)
        parsed (re-matches #".*epoch(\d+)-gen(\d+)\.pop" (.getName file))]
    (when parsed
      (let [[_ epoch gen] parsed]
        {:type ::fileinfo
         :file file,
         :epoch (Integer/parseInt epoch),
         :gen (Integer/parseInt gen)}))))

(defn saved-files [dir]
  (let [dir (mkfile dir)]
    (->> (.listFiles dir)
      (keep parse-filename))))

(defn epochs [dir]
  (->> dir saved-files (map :epoch) distinct sort))

(defn just-epoch [epoch file-infos]
  (->> file-infos
    (filter #(= epoch (:epoch %)))))

(defn just-gen [gen file-infos]
  (->> file-infos
    (filter #(= gen (:gen %)))
    first))

(defn just-nth [n file-info]
  (-> file-info :file edn-slurp (nth n)))

(defn zeroth-gen-of-epoch [epoch file-infos]
  (->> file-infos
    (some #(and (= epoch (:epoch %))  ;TODO Call just-epoch and just-gen
                (zero? (:gen %))))))

(defn last-gen-of-epoch [epoch file-infos]
  (->> file-infos
    (filter #(= epoch (:epoch %)))  ;TODO Call just-epoch
    (apply max-key :gen)))

(defn read-gt
 ([dir epoch gen]
  (read-gt dir epoch gen 0))  ;index 0 is the best genotype
 ([dir epoch gen index]
  (->> dir
    saved-files
    (just-epoch epoch)
    (just-gen gen)
    (just-nth index))))

(defn best-of-gen [file-info]
  (-> file-info :file readbest))

(defn save-vfn-of [epoch gen {:keys [data-dir] :as opts}]
  (let [gt (read-gt data-dir epoch gen)
        file (make-file (str "epoch" epoch "-gen" gen ".gfn") opts)]
    (save-fn (vfn gt) file)))

(defn save-gfn-of [epoch gen {:keys [data-dir] :as opts}]
  (let [gt (read-gt data-dir epoch gen)
        file (make-file (str "epoch" epoch "-gen" gen ".gfn") opts)]
    (save-fn (gfn gt) file)))

(defn save-wfn-of [epoch gen {:keys [data-dir] :as opts}]
  (let [gt (read-gt data-dir epoch gen)
        file (make-file (str "epoch" epoch "-gen" gen ".gfn") opts)]
    (save-fn (wfn gt) file)))

(defn pr-epoch-data [dir epoch]
  (let [best-gt (->> dir saved-files (last-gen-of-epoch epoch) best-of-gen)]
    (apply println (mround-floats (genotype-data best-gt)))))

(defn data-epoch-by-epoch [dir]
  (doseq [epoch (epochs dir)]
    (pr-epoch-data dir epoch)))

(defn pr-gen-data [dir epoch gen]
  (apply println gen (-> (read-gt dir epoch gen) genotype-data mround-floats)))

(defn data-gen-by-gen [dir epoch]
  (doseq [gen (->> dir saved-files (just-epoch epoch) (sort-by :gen))]
    (let [best-gt (best-of-gen gen)]
      (apply println (mround-floats (genotype-data best-gt))))))

(defn save-data-gen-by-gen [epoch {:keys [data-dir] :as opts}]
  (let [datafile (io/file data-dir (str "epoch" epoch ".data"))]
    (io/make-parents datafile)
    (with-open [datafile (io/writer datafile)]
      (with-*out* datafile
        (data-gen-by-gen data-dir epoch)))))

(defn save-conclusion-of-epoch [epoch {:keys [data-dir] :as opts}]
  (let [datafile (io/file data-dir (str "epochs.data"))]
    (io/make-parents datafile)
    (with-open [datafile (io/writer datafile :append true)]
      (with-*out* datafile
        (pr-epoch-data data-dir epoch)))))

;;; futures ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(def futures (atom {}))

(defn time-for-task? [{:keys [epoch generation] :as task} population opts]
  (let [task (if (= :last (:epoch task))
               (assoc task :epoch (:n-epochs opts))
               task)
        task (if (= :last (:generation task))
               (assoc task :generation (:generations opts))
               task)
        {:keys [epoch generation]} task]
    (cond
      (nil? epoch)
        (when (:all-epochs-done? population)
          task)
      (not= epoch (:epoch population))
        nil
      (nil? generation)
        (when (:epoch-done? population)
          task)
      (= generation (:generation population))
        task
      nil)))

(defn all-futures-done? [{:keys [data-dir] :as opts}]
  (let [is? #(= data-dir (:data-dir %))]
    (every? :done?
      (S/select [(S/selected? S/MAP-KEYS is?) S/MAP-VALS] @futures))))

(defn await-futures [opts]
  (when (not (all-futures-done? opts))
    (Thread/sleep 500)
    (recur opts)))

(def future-types
  {:save-vfn
   (fn [task-k population opts]
     (let [vfnfile (task-k->file task-k "-best.vfn" opts)
           gfnfile (task-k->file task-k "-best.gfn" opts)
           dotfile (task-k->file task-k "-best.dot" opts)]
       (save-dot dotfile (best-of population) opts)
       (save-vfn (best-of population) vfnfile)
       (matlab (str "mesh3 '" vfnfile "'"))
       (save-gfn (best-of population) gfnfile)
       (matlab (str "scat3 '" gfnfile "'"))))
   :epochs-data
     (fn [task-k population opts]
       (let [datafile (task-k->file task-k "epochs.data" opts)]
         (with-open [datafile (io/writer datafile)]
           (with-*out* datafile
             (data-epoch-by-epoch (:data-dir task-k))))
         (matlab (str "epochs '" datafile "'"))))
     })

(defn future-id [task]
  (select-keys task [:data-dir :epoch :generation :type]))

(defn mark-future-done [id]
  (swap! futures assoc-in [id :done?] true))

(defn start-future
  [{:keys [type] :as task-k} population {:keys [data-dir] :as opts}]
  (let [task-k (merge task-k {:data-dir data-dir})
        f (get future-types type)
        id (future-id task-k)]
    (swap! futures assoc id {:done? false
                             :future (future
                                       (let [result (f task-k population opts)]
                                         (mark-future-done id)
                                         result))})))

(defn start-future-if-time [task-k population opts]
  (when-let [task-k (time-for-task? task-k population opts)]
    (start-future task-k population opts)))
  
(defn start-scheduled-futures [population opts]
  (doseq [task-type (keys future-types)]
    (let [schedule (get opts task-type)]
      (cond
        (map? schedule)
          (start-future-if-time
            (assoc schedule :type task-type) population opts)
        (vector? schedule)
          (doseq [task-k schedule]
            (start-future-if-time
              (assoc task-k :type task-type) population opts))
        schedule
          (start-future-if-time {:type task-type} population opts)
        nil))))

;;; running ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn write-done [{:keys [data-dir] :as opts}]
  (let [donefile (make-file "DONE" opts)]
    (with-open [donefile (io/writer donefile)]
      :nothing)))

;(println (->> "seed1-d2" saved-files (last-gen-of-epoch 5) best-of-gen))

;TODO Run in a higher dimension.

(defn run-epoch [start-population epoch & opts]
  (let-ga-opts opts
    (let [[w w-source] (make-epoch-w)]
      (println w-source)
      (with-state [p start-population]
        (assoc :epoch epoch, :generation 0, :epoch-done? false,
               :w-source w-source, :w (memoize w))
        (when (empty? (:individuals p))
          (make-random-population population-size))
        (accumulate-data)
        -- (save-gen-data p)
        -- (start-scheduled-futures p opts)
        ;-- (println (best-fitness-of p))
        ;-- (apply println (-> p best-of genotype-data mround-floats))
        (doseq [generation (range 1 (inc generations))]
          (assoc :generation generation)
          (next-gen)
          (accumulate-data)
          -- (save-gen-data p)
          ;-- (println (best-fitness-of p))
          ;-- (apply println (-> p best-of genotype-data mround-floats))
          -- (start-scheduled-futures p opts)
          )
        ;-- (save-dot (str "epoch" epoch "-best.dot") (best-of p) opts)
        (assoc :epoch-done? true)
        -- (start-scheduled-futures p opts)
;        (when (:save-epoch-vfn opts)
;          -- (start-future (save-vfn (best-of p)
;                                     (make-file (str "epoch" epoch "-best.vfn")
;                                                opts))))
;        -- (start-future (save-data-gen-by-gen epoch opts))
;        -- (start-future (save-conclusion-of-epoch epoch opts))
        ))))

(defn run [& opts]
  (let-ga-opts opts
    (binding [*data-dir* data-dir
              *genotype->phenotype* genotype->phenotype]
      (with-rng-seed seed
        (println (str "seed=" seed))
        (clear-data-dir data-dir)
        (save-opts opts)
        (with-state [p empty-population]
          (assoc :seed util/*rng-seed*, :all-epochs-done? false,
                 :epoch 0, :generation 0)
          -- (start-scheduled-futures p opts)
          (doseq [epoch (range 1 (inc n-epochs))]
            (when (> n-epochs 1)
              -- (print (str "epoch" epoch \space)))
            (run-epoch epoch opts)
            ;-- (apply println (-> p best-of genotype-data mround-floats)))
            )
          (assoc :all-epochs-done? true)
          ;-- (save-dot "best.dot" (best-of p) opts)
          -- (start-scheduled-futures p opts)
          -- (await-futures opts)
          -- (write-done opts)
          (return (best-of p)))))))

#_(defn run-and-save [& opts]
  (let [opts (->> opts (map clojure.edn/read-string) (apply hash-map))
        winning-genotype (run opts)]
    (spit "winner" winning-genotype)
    (spit "winner.dot" (dot winning-genotype))
    (with-*out* (io/writer "winner.acclivity")
      (let [climber-results (hill/run-climbers (vfn winning-genotype))]
        (run! println climber-results)
        (println (util/average (map first climber-results)))))
    (save-vfn winning-genotype "winner.vfn")))

(defn test-run [& opts]
  (let [opts {:seed 0 :population-size 40 :generations 20 :n-epochs 100
              :save-vfn [{:epoch 1 :generation 0}
                         {:epoch :last :generation :last}]
              :epochs-data true}]
    (time
      (doseq [seed [1 2 3 4 5]]
        (run (assoc opts :seed seed))))))

#_(defn scripted-run []
  (run :seed 1 :population-size 40 :generations 20 :n-epochs 40
       :save-vfn [{:epoch 1 :generation 0} {:epoch 20 :generation :last}
                  {:epoch 20 :generation :last}]
       :save-gfn [{:epoch 1 :generation 0} {:epoch 20 :generation :last}
                  {:epoch 20 :generation :last}]))

(defn make-article-data
  "Makes all the fitness-func and population data files for the article."
  []
  (doseq [seed (range 1 41)]
    (run :seed seed)))

(defn -main [& args]
  (println "For now, run in the REPL."))
