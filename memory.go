package sparse

import "github.com/Konstantin8105/mms"

// TODO: check data race
var (
	floats = mms.Float64sCache{}
	ints   = mms.IntsCache{}
)

// func init() {
// 	mms.Debug = true
// }
