(ns farg.acclivation
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.math.combinatorics :as combo]
            [clojure.math.numeric-tower :as math]
            [clojure.java.io :refer [file writer]]
            [clj-time.local :as ltime]
            [farg.util :refer [dd dde choose choose-one with-rng-seed mround
                               defopts with-*out*]
             :as util]
            [farg.with-state :refer [with-state]]
            [farg.acclivation.sa :as sa]
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
; run-epoch

(def ^:dynamic *w*)  ;phenotype fitness
(def ^:dynamic *genotype->phenotype*)
(def ^:dynamic *pop-file* *out*) ;output file for population

(defn strip-type [x]
  (vary-meta x dissoc :type))

(defn graph-cstr [graph]
  (str (sort (uber/nodes graph)) \space
       (uber/count-edges graph) " edges"))

(defn genotype-cstr [g]
  (str (-> (strip-type g)
           (select-keys [:numbers :phenotype :fitness])) \space
       (graph-cstr (:graph g))))

;  (let [graph (:graph g)
;        g-numbers (:numbers g)
;        ph-numbers (:phenotype g)]
;    (str (:numbers g) \space
;         (util/nilstr ph-numbers) \space
;         (graph-cstr graph))))

(defn population-cstr [p]
  (str (strip-type (select-keys p [:generation :epoch])) \newline
       (count (:individuals p))))

(defn convenient-str ^String [x]
  (if (map? x)
      (case (:type x)
        :genotype (genotype-cstr x)
        :population (population-cstr x)
        (apply str (strip-type x)))
      (apply str (strip-type x))))

(defmethod print-method :acclivation [v ^java.io.Writer w]
  (.write w (convenient-str v)))

(defn add-phenotype [g]
  (if (contains? g :phenotype)
      g
      (assoc g :phenotype (*genotype->phenotype* g))))

(defn add-fitness [g]
  (let [g (add-phenotype g)]
    (if (contains? g :fitness)
        g
        (assoc g :fitness (*w* (:phenotype g))))))

(defn genotype-str [g]
  (->> (add-fitness g)
       (genotype-cstr)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

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
  (let [invv (make-inverted-v 1.50 0.5)]
    (fn [[x1 x2]]
      (if (zero? x1)
          0.0
          (invv (/ x2 x1))))))

(def w2
  (let [near1 (util/piecewise-linear 0.0 0.0, 1.0 1.0, 1.0 0.0)]
    (fn [[x1 x2]]
      (* 5 (+ (near1 x1) (near1 x2))))))

(defn w12 [ph]
  (* (w-ratio ph) (w2 ph)))

(defn w-many-small-hills [[x1 x2]]
  (* 2 (Math/cos (* 30 x1))
       (Math/sin (* 30 x2))))

(defn w-distance [[x1 x2]]
  (let [center [0.2 0.2]
        invv (make-inverted-v 0.0 0.5)]
    #_(* 5.0 (invv (util/distance [x1 x2] center)))
    (- 2.0 (util/distance [x1 x2] center))
    ))

(defn w [ph]
  (+ (w12 ph) (w-many-small-hills ph) (w-distance ph)))

;(def w (memoize w))

#_(defn w [ph]
  (* #_(w-equal ph) (Math/pow (w-distance ph) 2.0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-node [n]
  {:pre [(integer? n)]}
  (keyword (str \n n)))

(defn random-weight []
  (choose-one -1.0 1.0))

(defn make-random-graph []
  (let [n-extra-nodes (util/rand-int 2 10)
        extra-nodes (->> (take n-extra-nodes (iterate inc 1))
                         (map make-node))
        nodes (into [:g1 :g2 :p1 :p2] extra-nodes)
        n-edges (util/rand-int 5 (* 3 (count nodes)))]
    (with-state [g (apply uber/digraph nodes)]
      (dotimes [n n-edges]
        (bind [from to] (util/choose-with-replacement 2 nodes))
        (bind weight (random-weight))
        (uber/add-edges [from to weight]))
      ;(uber/add-edges [:g1 (choose extra-nodes) (random-weight)])
      ;(uber/add-edges [:g2 (choose extra-nodes) (random-weight)])
      ;(uber/add-edges [(choose extra-nodes) :p1 (random-weight)])
      ;(uber/add-edges [(choose extra-nodes) :p2 (random-weight)])
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
   :epoch 1})

(defn make-random-population [n]
  (assoc empty-population
         :individuals (vec (take n (repeatedly make-random-genotype)))
         :population-size n))

(defn extract-phenotype [g]
  [(sa/a g :p1) (sa/a g :p2)])

(defn genotype->phenotype [g]
  (-> (:graph g)
      (sa/spread-activation (zipmap [:g1 :g2] (:numbers g)) :iterations 20)
      (extract-phenotype)))

;(def genotype->phenotype (memoize genotype->phenotype))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn turn-knob [g]
  (let [i (choose [0 1])
        Δ (util/rand -0.02 0.02)]
    (update-in g [:numbers i] #(sa/squash (+ % Δ)))))

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
  (update g :graph #(uber/remove-edges % (choose (uber/edges %)))))

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
  [[turn-knob 8]
   [move-edge 4]
   [add-node 1]
   [remove-node 1]
   [add-edge 1]
   [remove-edge 1]])

(defn mutate [g]
  ((util/weighted-choice mutations) g))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn mutate-n [n g]
  (take n (repeatedly #(mutate g))))

(defn mutate-population [n p]
  (assoc p :individuals (->> (:individuals p)
                             (mapcat #(mutate-n n %)))))

(defn take-first-half [coll]
  (let [coll (seq coll)]
    (take (/ (count coll) 2) coll)))

(defn take-second-half [coll]
  (let [coll (seq coll)]
    (drop (/ (count coll) 2) coll)))

(defn crossover [g1 g2]
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

(defn random-pair [coll]
  (util/choose-without-replacement 2 coll))

(defn mutate-and-crossover [mutate crossover n-by-mutation n-by-crossover p]
  (let [parents (:individuals p)
        mutants (->> (repeatedly #(choose parents))
                     (mapcat #(mutate %))
                     (take n-by-mutation))
        crosses (->> (repeatedly #(random-pair parents))
                     (mapcat #(apply crossover %))
                     (take n-by-crossover))]
    (assoc p :individuals (into mutants crosses))))

(defn tournament-selection
  [f-fitness n p]
  (apply max-key f-fitness
         (util/choose-without-replacement n (:individuals p))))
;  [f-fitness n p]
;  (let [indices (util/choose-without-replacement n
;                  (range (count (:individuals p))))
;        individuals (add-fitnesses f-fitness (:individuals p) indices)
;        contestants (map individuals indices)]
;    [(assoc p :individuals individuals)
;     (apply max-key :fitness contestants)]))

(defn tournament-selection-on-population [f-fitness n p]
  (loop [winners #{}]
    (if (>= (count winners) (:population-size p))
        (assoc p :individuals winners)
        (recur (conj winners (tournament-selection f-fitness n p))))))

(defn next-generation [vary select population]
  (-> population
      vary
      select
      (update :generation inc)))

(defn genotype-fitness [g]
  (->> (genotype->phenotype g)
       (w)))

(defn default-to [x default]
  (if (some? x) x default))

(defn dot-edge [graph edge]
  (str (name (:src edge)) " -> " (name (:dest edge)) \space
       "[label=" (mround (uber/weight graph edge)) "]; "))
       ;"[label=" (mround (sa/a graph edge)) "]; "))

(defn dot-node [graph node]
  (str (name node) \space "[label=" (mround (sa/a graph node)) "] "))

(defn make-dot [graph]
  (str "digraph graphname {
{ rank=source edge [style=\"invis\"] g1 -> g2 }
{ rank=sink edge [style=\"invis\"] p1 -> p2 }
"
       (apply str (map #(dot-node graph %) (uber/nodes graph)))
       (apply str (map #(dot-edge graph %) (uber/edges graph)))
       "}"))

(defn print-genotype-graph [g]
  (let [graph (sa/spread-activation (:graph g)
                                    (zipmap [:g1 :g2] (:numbers g))
                                    :iterations 20)]
    (println (make-dot graph))))

(defn print-population [p]
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn make-fitness-data []
  (doseq [x (range -1.01 1.01 0.005)
          y (range -1.01 1.01 0.005)]
    (println x y (w [x y]))))

(defn trim-round-off-error [x]
  (/ (Math/round (* 100000000.0 x))
     100000000.0))

(defn fitness-as-seen-by [g & {:keys [radius step]
                               :or {radius nil, step 0.005}}]
  ;(let [[g1 g2] [-1.01 1.01]]  ;(:numbers g)]  ;or else [-1.01 1.01]
  (let [[g1 g2] (:numbers g)
        xrange (if radius (range (- g1 radius) (+ g1 radius) step)
                          (range -1.0 1.000000000001 step))
        yrange (if radius (range (- g2 radius) (+ g2 radius) step)
                          (range -1.0 1.000000000001 step))]
    (doseq [x xrange, y yrange]
      (let [x (trim-round-off-error x)
            y (trim-round-off-error y)
            g' (assoc g :numbers [x y])
            ph (genotype->phenotype g')
            [phx phy] ph
            fitness (w ph)]
        (println x y fitness phx phy)))))

(defn print-fitness-fn
  "(f [startcoord]) should return [endcoord fitness]."
  [f & {:keys [step] :or {step 0.005}}]
  (let [xrange (range -1.0 1.000000000001 step)
        yrange (range -1.0 1.000000000001 step)]
    (doseq [x xrange, y yrange]
      (let [startx (trim-round-off-error x)
            starty (trim-round-off-error y)
            startcoord [startx starty]
            [endcoord fitness] (f [startx starty])]
        (apply println (concat startcoord [fitness] endcoord))))))

(defn virtual-fitness-fn
  "Returns the fitness function f seen by the \"numbers\" part of genotype g.
  (f [startcoord]) returns [endcoord fitness], where endcoord is the phenotype
  resulting from startcoord."
  [g]
  (fn [startcoord]
    (let [g' (assoc g :numbers startcoord)
          ph (genotype->phenotype g')]
      [ph (w ph)])))

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

(defn print-w
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn now-suffix []
  (ltime/format-local-time (ltime/local-now) :basic-date-time-no-ms))

(defopts let-ga-opts
  {:keys [generations population-size n-mutants n-crosses f-mutate
                     f-cross tourney-size vary select fitness seed
                     radius step]
   :or {generations 20, population-size 20, tourney-size 5,
        f-mutate (partial mutate-n 2), f-cross crossover,
        fitness genotype-fitness, seed 1,
        radius nil, step 0.01}} ;arguments for fitness-as-seen-by
          ;make step 0.005 for precise fitness func (it just takes a long time)
  n-mutants (default-to n-mutants (int (* population-size 1.5)))
  n-crosses (default-to n-crosses (int (* population-size 0.5)))
  vary (default-to vary (partial mutate-and-crossover f-mutate f-cross
                                 n-mutants n-crosses))
  select (partial tournament-selection-on-population fitness tourney-size)
  next-gen (partial next-generation vary select)
  prefix (str "seed" seed "-") ;prefix for output filenames
  pop-file (file (str prefix "pop"))
    ;population of each generation
  vfit-file (file (str prefix "vfit"))
    ;fitness-as-seen-by best of last gen
  )

(defn print-phenotype-fitness-fn [& opts]
  (let-ga-opts opts
    (with-*out* (writer "phenotype-fitness")
      (print-fitness-fn (fn [ph] [ph (w ph)])))))

;(defn run [& {:keys [generations population-size n-mutants n-crosses f-mutate
;                     f-cross tourney-size vary select fitness seed]
;              :or {generations 20, population-size 40, tourney-size 5,
;                   f-mutate (partial mutate-n 2), f-cross crossover,
;                   fitness genotype-fitness, seed 1}}]
;  (let [n-mutants (default-to n-mutants (int (* population-size 1.5)))
;        n-crosses (default-to n-crosses (int (* population-size 0.5)))
;        vary (default-to vary (partial mutate-and-crossover f-mutate f-cross
;                                       n-mutants n-crosses))
;        select (partial tournament-selection-on-population fitness tourney-size)
;        next-gen (partial next-generation vary select)]
(defn run [& opts]
  (let-ga-opts opts
    (println "Sending output to" (.getName pop-file) (.getName vfit-file) "...")
    (binding [*w* w
              *genotype->phenotype* genotype->phenotype
              *pop-file* (writer pop-file)]
      (with-rng-seed seed
        (with-state [p (make-random-population population-size)]
          -- (print-population p)
          (dotimes [i generations]
            (next-gen)
            -- (print-population p))
          (bind best (best-of p))
          -- (with-*out* (writer vfit-file)
               (print-fitness-fn (virtual-fitness-fn best) :step step))
;          -- (with-*out* (writer vfit-file)
;               (fitness-as-seen-by best :radius radius :step step))
          (return best))))))

(defn make-article-data
  "Makes all the fitness-func and population data files for the article."
  []
  (doseq [seed (range 1 41)]
    (run :seed seed)))

(defn -main [& args]
  (println "For now, run in the REPL."))
