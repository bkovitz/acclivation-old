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
            [farg.util :as util :refer [dd dde]]
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

(defn hill-step
  "Returns point with highest f within one step of xx. If none have higher
  f than xx, returns xx."
  [f step xx]
  (let [new-xx (apply max-key f (all-steps step xx))]
    (if (> (f new-xx) (f xx)) new-xx xx)))

(defn hill-climb
  "Hill-climbs 'f', one 'step' at a time, starting at 'xx'. Returns
  [best-fitness start-xx best-xx]."
  [f step start-xx]
  (loop [best-fitness (f start-xx), xx start-xx]
    (let [new-xx (hill-step f step xx)]
      (if (= xx new-xx)
        [(normalize best-fitness) start-xx xx]
        (recur (f new-xx) new-xx)))))

(defn run-climbers
  [f step dimension n-climbers]
  (->> (deterministic-random-xxs step dimension)
       (take n-climbers)
       (map #(hill-climb f step %))))
