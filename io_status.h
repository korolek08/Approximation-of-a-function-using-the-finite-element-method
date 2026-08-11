typedef enum io_status_ {
	SUCCESS,
	ERROR_OPEN,
	ERROR_READ,
	ERROR_FORMAT,
	ERROR_UNKNOWN,
	DEGENERATE,
	PROCESSING,
	READY,
	START,
	FINISH,
	undef,
} io_status;

enum task_status {
  WAIT,
  DO,
  STOP,
};

enum calculation {
  PROCESS,
  END,
  REST,
};
