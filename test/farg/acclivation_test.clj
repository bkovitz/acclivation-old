(ns farg.acclivation-test
  (:refer-clojure :exclude [rand rand-int cond memoize])
  (:require [better-cond.core :refer [cond]]
            [clojure.edn :as edn]
            [clojure.tools.trace :refer [deftrace] :as trace]
            [clojure.pprint :refer [pprint]]
            [clojure.test :refer :all]
            [farg.acclivation :refer :all]
            [farg.util :as util :refer [dd dde]]
            [farg.with-state :refer [with-state]]))

(deftest test-edn
  (let [gt (make-random-genotype)
        sgt (genotype->edn gt)]
    (let [ht (edn/read-string {:readers edn-readers} sgt)]
      (is (= gt ht)))))

(deftest test-run
  (run :seed 0 :population-size 5 :generations 2 :n-epochs 2))
