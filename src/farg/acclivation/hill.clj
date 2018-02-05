(ns farg.acclivation.hill
  "Hill-climbing, to measure acclivity of a function.

  All coordinates are in the range [-1.0, +1.0].

  Naming convention: xx is a vector of coordinates. x is a single number."
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.math.combinatorics :as combo]
            [clojure.math.numeric-tower :as math]
            [clojure.java.io :refer [file writer]]
            [clj-time.local :as ltime]
            [farg.util :as util :refer [dd dde with-rng-seed]]
            [farg.with-state :refer [with-state]]))

(defn normalize [x]
  (if (zero? x)
    0.0  ;force 0.0; prevent -0.0
    (-> x (* 100000) (Math/round) (/ 100000.0))))

(defn lazy-deterministic-shuffled-indices [seed ub]
  (let [rng (java.util.Random. seed)]
    (letfn [(next-int [] (.nextInt rng ub))
            (lazy [so-far]
              (lazy-seq
                (when (< (count so-far) ub)
                  (let [i (->> (repeatedly next-int)
                               (filter #(not (contains? so-far %)))
                               (first))]
                    (cons i (lazy (conj so-far i)))))))]
      (lazy #{}))))

(defn make-flatten-coords
  "Returns a function that maps an integer index to a multidimensional
  coordinate.

  row: A vector of all the numbers in a single dimension.

  dimension: The number of dimensions in the coordinates."
  [row dimension]
  (let [row-size (count row)]
    (fn [index]
      (loop [index index, coords [], n dimension]
        (cond
          (zero? n) coords
          :let [i (rem index row-size)
                index (quot index row-size)]
          (recur index (conj coords (get row i)) (dec n)))))))

(defn steps-1d [step]
  (->> (range -1.0 1.0 step)
    (map normalize)
    vec))

(defn deterministic-random-xxs [step dimension]
  (let [row (steps-1d step)
        i->xx (make-flatten-coords row dimension)
        n-indices (int (Math/pow (count row) dimension))]
    (->> (lazy-deterministic-shuffled-indices 1 n-indices)
         (map i->xx))))

(defn in-range? [xx]
  (every? #(<= -1.0 % 1.0) xx))

(defn all-steps
  "Returns all coordinates that vary from 'xx' by 'step' in one dimension."
  [step xx]
  (->> (for [n (range (count xx))
             direction [- +]]
         (update xx n #(normalize (direction % step))))
       (filter in-range?)))

(defn +vector [xx yy]
  (->> (map #(normalize (+ %1 %2)) xx yy)
       vec))

#_(defn all-steps
  "all-steps for steepest ascent (diagonal moves allowed)."
  [step xx]
  (->> (repeat (count xx) [(- step) 0.0 (+ step)])
       (apply combo/cartesian-product)
       (map #(+vector % xx))
       (remove #(every? zero? %))
       (filter in-range?)))

(defn- update-best [f [bests best-fitness :as prev] xx]
  (cond
    :let [fitness (f xx)]
    (empty? bests)
      [#{xx} fitness]
    (< fitness best-fitness)
      prev
    (= fitness best-fitness)
      [(conj bests xx) best-fitness]
    [#{xx} fitness]))

(defn hill-step
  "Returns point with highest f within one step of xx. If there's a tie
  between xx and another point, returns the other point (chosen arbitrarily
  but deterministically if more than one other point has the same fitness)."
  [f step xx]
  (with-rng-seed 1
    (let [start-fitness (f xx)]
      (letfn [(update-best [[bests best-fitness :as prev] xx]
                (cond
                  :let [fitness (f xx)]
                  (< fitness start-fitness)
                    prev
                  (empty? bests)
                    [#{xx} fitness]
                  (< fitness best-fitness)
                    prev
                  (= fitness best-fitness)
                    [(conj bests xx) best-fitness]
                  [#{xx} fitness]))]
        (let [[better-neighbors _] (reduce update-best
                                           [nil nil]
                                           (all-steps step xx))]
          (if (empty? better-neighbors)
            xx  ;local optimum
            (util/choose-from better-neighbors)))))))

(defn hill-climb
  "Hill-climbs 'f', one 'step' at a time, starting at 'xx'. Stops climbing
  when reaching a local optimum, or, on a neutral plateau, when either
  the same point has been reached 5 times or fitness has not improved for
  1000 steps.  Returns [best-fitness n-steps start-fitness start-xx best-xx]."
  [f step start-xx]
  (let [start-fitness (f start-xx)]
    (loop [best-fitness start-fitness
           n-steps 0
           xx start-xx
           seen-before {}
           n-same-fitness 0]
      (let [new-xx (hill-step f step xx)]
        (cond
          (= xx new-xx) ;hill-step didn't move: local optimum
            [best-fitness n-steps start-fitness start-xx xx]
          :let [new-fitness (f new-xx)]
          (not= new-fitness best-fitness)
            (recur new-fitness (inc n-steps) new-xx {} 0)
          (or (>= (get seen-before new-xx 0) 100)
              (>= n-same-fitness 1000))
            [best-fitness n-steps start-fitness start-xx new-xx]
          (recur best-fitness
                 (inc n-steps)
                 new-xx
                 (update seen-before new-xx (fnil inc 0))
                 (inc n-same-fitness)))))))

(defn run-climbers
  [f step dimension n-climbers]
  (->> (deterministic-random-xxs step dimension)
       (take n-climbers)
       (map #(hill-climb f step %))))
