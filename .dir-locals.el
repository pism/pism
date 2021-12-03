;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

((nil
  (fill-column . 90)
  (eval . (condition-case nil
              (progn
                (add-to-list 'grep-find-ignored-files "*TAGS" )
                (ws-butler-mode))
            (error nil)))
  (eval . (progn
            (c-set-offset 'innamespace 0)
            (c-set-offset 'inlambda 0)))))
