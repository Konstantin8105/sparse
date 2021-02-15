package sparse

import (
	"fmt"
	"strings"
)

// Errors of work
var (
	ErrInputData      = wrapError{msg: "input data"}
	ErrInitialization = wrapError{msg: "initialization"}
	ErrAllocation     = wrapError{msg: "allocation"}
	ErrReallocation   = wrapError{msg: "realloc"}
	ErrCriteria       = wrapError{msg: "criteria"}
)

type wrapError struct {
	err error
	msg string
}

func (err wrapError) Error() string {
	if err.err != nil {
		return fmt.Sprintf("%s: %v", err.msg, err.err)
	}
	return err.msg
}

func (err wrapError) wrap(inner error) error {
	return wrapError{msg: err.msg, err: inner}
}

func (err wrapError) Unwrap() error {
	return err.err
}

func (err wrapError) Is(target error) bool {
	ts := target.Error()
	return ts == err.msg || strings.HasPrefix(ts, err.msg+": ")
}
