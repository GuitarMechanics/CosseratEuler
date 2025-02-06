import cmd

class NumberGame(cmd.Cmd):
    intro = '숫자게임에 오신것을 환영합니다. 도움말은 help 또는 ? 을 입력하세요.\n'
    prompt = '(숫자게임) '

    def do_quit(self, *arg):
        return True


if __name__ == '__main__':
    NumberGame().cmdloop()